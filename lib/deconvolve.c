/*********************************************************************
Functions for deconvolution operations.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Alvaro Morales <alvaro96m@hotmail.com>
Contributing author(s):
Copyright (C) 2024-2024 Free Software Foundation, Inc.

Gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <errno.h>
#include <error.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <string.h>

#include <gnuastro/complex.h>
#include <gnuastro/deconvolve.h>
#include <gnuastro/fft.h>
#include <gnuastro/pointer.h>

/**
 * @brief Implement a deconvolution using the Wiener / Tikhonov method
 * described at https://ui.adsabs.harvard.edu/abs/2002PASP..114.1051S/abstract
 * (eq 7)
 *
 * @param image The image with distorsion and noise
 * @param PSF The distorsion kernel
 * @param lambda
 * @param numthreads Number of threads in FFT
 * @param output
 */
void
gal_deconvolve_tikhonov (const gal_data_t *image, const gal_data_t *PSF,
                         double lambda, size_t numthreads, gal_data_t **output)
{
  gsl_complex_packed_array imagepadding;  // original image after padding
  gsl_complex_packed_array psfpadding;    // kernel after padding
  gsl_complex_packed_array psffreq;       // PSF(u,v)
  gsl_complex_packed_array imagefreq;     // I(u,v)
  gsl_complex_packed_array psffconj;      // PSF*(u,v)
  gsl_complex_packed_array psffreqsquare; // |PSF(u,v)|^2
  gsl_complex_packed_array numerator;     // PSF*(u,v)I(u,v)
  gsl_complex_packed_array denominator;   // |PSF(u,v)|^2 + λ
  gsl_complex_packed_array
      deconvolutionfreq; // PSF*(u,v)I(u,v)/ (|PSF(u,v)|^2 + λ)
  gsl_complex_packed_array deconvolution; // Deconvolution in time domine

  gal_data_t *data = NULL;
  size_t dsize[2];
  size_t size;
  double *tmp;

  /* Check image type. */
  if (image->type != GAL_TYPE_FLOAT32)
    error (EXIT_FAILURE, 0, "%s: input data must be float 64", __func__);

  /* Process the image and kernel to have same size and be complex numbers. */
  gal_complex_create_padding (image, PSF, &imagepadding, &psfpadding, dsize,
                              &dsize[1]);
  size = dsize[0] * dsize[1]; // Total number of elements.

  /* Normalize and rearange the kernel. */
  gal_complex_normalize (psfpadding, size);
  gal_fft_swap_quadrant (psfpadding, dsize);

  /* Convert to frequency domain. */
  gal_fft_two_dimension_transformation (psfpadding, dsize, &psffreq,
                                        numthreads, gsl_fft_forward);
  gal_fft_two_dimension_transformation (imagepadding, dsize, &imagefreq,
                                        numthreads, gsl_fft_forward);

  /* Calculate numerator PSF*(u,v)I(u,v) */
  gal_complex_conjugate (psffreq, size, &psffconj);
  gal_complex_multiply (psffreq, imagefreq, &numerator, size);

  /* Caculate denominator |PSF(u,v)|^2 + λ */
  gal_complex_multiply (psffreq, psffconj, &psffreqsquare, size);
  gal_complex_add_scalar (psffreqsquare, size, lambda + I * 0, &denominator);

  /* Calculate the deconvolve image (in frequency domain).*/
  gal_complex_divide (numerator, denominator, &deconvolutionfreq, size,
                      lambda);

  /* Go back to time domain. */
  gal_fft_two_dimension_transformation (
      deconvolutionfreq, dsize, &deconvolution, numthreads, gsl_fft_backward);

  /* Convert to Real number and convert it to GAL TYPE.*/
  gal_complex_to_real (deconvolution, dsize[0] * dsize[1],
                       COMPLEX_TO_REAL_REAL, &tmp);
  data
      = gal_data_alloc (tmp, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1, MIN_MAP_SIZE,
                        1, NULL, NULL, NULL); // data has to be 32
  *output = data;

  /* Free resources. */
  free (imagepadding);
  free (psfpadding);
  free (psffreq);
  free (psffconj);
  free (imagefreq);
  free (psffreqsquare);
  free (numerator);
  free (denominator);
  free (deconvolutionfreq);
  free (deconvolution);
}
