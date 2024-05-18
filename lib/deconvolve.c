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
                         double lambda, size_t numthreads, size_t minmapsize,
                         gal_data_t **output)
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
    error (EXIT_FAILURE, 0, "%s: input data must be float 32", __func__);

  /* Process the image and kernel to have same size and be complex numbers. */
  gal_complex_create_padding (image, PSF, &imagepadding, &psfpadding, dsize,
                              &dsize[1]);
  size = dsize[0] * dsize[1]; // Total number of elements.

  /* Normalize and rearange the kernel. */
  gal_complex_normalize (psfpadding, size);
  gal_fft_shift_center (psfpadding, dsize);

  /* Convert to frequency domain. */
  gal_fft_two_dimension_transformation (
      psfpadding, dsize, &psffreq, numthreads, minmapsize, gsl_fft_forward);
  gal_fft_two_dimension_transformation (imagepadding, dsize, &imagefreq,
                                        numthreads, minmapsize,
                                        gsl_fft_forward);

  /* Calculate numerator PSF*(u,v)I(u,v) */
  psffconj = gal_complex_conjugate (psffreq, size);
  numerator = gal_complex_multiply (psffconj, imagefreq, size);

  /* Caculate denominator |PSF(u,v)|^2 + λ */
  psffreqsquare = gal_complex_multiply (psffreq, psffconj, size);
  denominator = gal_complex_add_scalar (psffreqsquare, size, lambda + I * 0);

  /* Calculate the deconvolve image (in frequency domain).*/
  deconvolutionfreq
      = gal_complex_divide (numerator, denominator, size, lambda);

  /* Go back to time domain. */
  gal_fft_two_dimension_transformation (deconvolutionfreq, dsize,
                                        &deconvolution, numthreads, minmapsize,
                                        gsl_fft_backward);

  /* Convert to Real number and convert it to GAL TYPE.*/
  tmp = gal_complex_to_real (deconvolution, dsize[0] * dsize[1],
                             COMPLEX_TO_REAL_REAL);
  data = gal_data_alloc (tmp, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1, minmapsize,
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

/**
 * @brief Implements a deconvolution using the naive method:
 * O(u,v) = I(u,v)/PSF(u,v)
 *
 * @param image
 * @param PSF
 * @param numthreads
 * @param minmapsize
 * @param output
 */
void
gal_deconvolve_naive (const gal_data_t *image, const gal_data_t *PSF,
                      size_t numthreads, size_t minmapsize,
                      gal_data_t **output)
{
  gsl_complex_packed_array imagepadding;      // original image after padding
  gsl_complex_packed_array psfpadding;        // kernel after padding
  gsl_complex_packed_array psffreq;           // PSF(u,v)
  gsl_complex_packed_array imagefreq;         // I(u,v)
  gsl_complex_packed_array deconvolutionfreq; // I(u,v) / PSF(u,v)
  gsl_complex_packed_array deconvolution;     // Deconvolution in time domine

  gal_data_t *data = NULL;
  size_t dsize[2];
  size_t size;
  double *tmp;

  /* Check image type. */
  if (image->type != GAL_TYPE_FLOAT32)
    error (EXIT_FAILURE, 0, "%s: input data must be float 32", __func__);

  /* Process the image and kernel to have same size and be complex numbers. */
  gal_complex_create_padding (image, PSF, &imagepadding, &psfpadding, dsize,
                              &dsize[1]);
  size = dsize[0] * dsize[1]; // Total number of elements.

  /* Normalize and rearange the kernel. */
  gal_complex_normalize (psfpadding, size);
  gal_fft_shift_center (psfpadding, dsize);

  /* Convert to frequency domain. */
  gal_fft_two_dimension_transformation (
      psfpadding, dsize, &psffreq, numthreads, minmapsize, gsl_fft_forward);
  gal_fft_two_dimension_transformation (imagepadding, dsize, &imagefreq,
                                        numthreads, minmapsize,
                                        gsl_fft_forward);

  /* Calculate the deconvolve image (in frequency domain).*/
  deconvolutionfreq = gal_complex_divide (imagefreq, psffreq, size, 1e-6);

  /* Go back to time domain. */
  gal_fft_two_dimension_transformation (deconvolutionfreq, dsize,
                                        &deconvolution, numthreads, minmapsize,
                                        gsl_fft_backward);

  /* Convert to Real number and convert it to GAL TYPE.*/
  tmp = gal_complex_to_real (deconvolution, dsize[0] * dsize[1],
                             COMPLEX_TO_REAL_REAL);
  data = gal_data_alloc (tmp, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1, minmapsize,
                         1, NULL, NULL, NULL); // data has to be 32
  *output = data;

  /* Free resources. */
  free (imagepadding);
  free (psfpadding);
  free (psffreq);
  free (deconvolutionfreq);
  free (deconvolution);
}

void
richardson_lucy_init_solution (gsl_complex_packed_array *output, size_t size)
{
  gsl_complex_packed_array out;

  /* Allocate the space for the output array. */
  out = gal_pointer_allocate (GAL_TYPE_COMPLEX64, size, 1, __func__,
                              "rlobject");

  for (size_t index = 0; index < size; index++)
    {
      out[index * 2] = 1;
    }
  *output = out;
}

void
richardson_lucy_calculate_next_solution (gsl_complex_packed_array *solution,
                                         gsl_complex_packed_array h,
                                         gsl_complex_packed_array h_f,
                                         gsl_complex_packed_array image,
                                         double alpha, size_t *dsize,
                                         size_t minmapsize, size_t numthreads)
{
  gsl_complex_packed_array solution_f;    // FFT(x)
  gsl_complex_packed_array yest_f;        // h×FFT(x(k))
  gsl_complex_packed_array yest;          // y_est := FFT−1(h×FFT(x))
  gsl_complex_packed_array division;      // y/y_est
  gsl_complex_packed_array divisionfreq;  // FFT(y/y_est)
  gsl_complex_packed_array bracket_f;     // h*xFFT(y/y_est)
  gsl_complex_packed_array bracket;       // FFT^-1(bracket_f)
  gsl_complex_packed_array bracket_alpha; // FFT^-1(bracket_f)^alpha
  gsl_complex_packed_array next_solution; // FFT^-1(bracket_f)^alpha·x
  size_t size = dsize[0] * dsize[1];

  /* Convert solution to frequency domain. */
  gal_fft_two_dimension_transformation (
      *solution, dsize, &solution_f, numthreads, minmapsize, gsl_fft_forward);

  yest_f = gal_complex_multiply (solution_f, h, size);

  /*Calculate y_est*/
  gal_fft_two_dimension_transformation (yest_f, dsize, &yest, numthreads,
                                        minmapsize, gsl_fft_backward);

  /*Calculate y/y_est and its FFT*/
  division = gal_complex_divide (image, yest, size,
                                 1e-6); // check min value issues

  gal_fft_two_dimension_transformation (
      division, dsize, &divisionfreq, numthreads, minmapsize, gsl_fft_forward);

  /* Calculate bracket value in Freq and Space*/
  bracket_f = gal_complex_multiply (divisionfreq, h_f, size);

  gal_fft_two_dimension_transformation (bracket_f, dsize, &bracket, numthreads,
                                        minmapsize, gsl_fft_backward);
  // alpha value
  bracket_alpha = gal_complex_power (bracket, alpha, size);
  /* Calculate next object */
  next_solution = gal_complex_multiply (bracket_alpha, *solution, size);

  /* Free all internal variables */
  free (*solution);
  free (solution_f);
  free (yest_f);
  free (yest);
  free (division);
  free (divisionfreq);
  free (bracket_f);
  free (bracket);
  free (bracket_alpha);

  *solution = next_solution;
}

void
gal_deconvolve_richardson_lucy (const gal_data_t *image, const gal_data_t *PSF,
                                size_t iterations, double alpha,
                                size_t minmapsize, size_t numthreads,
                                gal_data_t **output)
{
  gsl_complex_packed_array imagepadding; // original image after padding I(x,y)
  gsl_complex_packed_array psfpadding;   // kernel after padding PSF(x,y)
  gsl_complex_packed_array psffreq;      // PSF(u,v)
  gsl_complex_packed_array psffconj;     // PSF*(u,v)
  gsl_complex_packed_array object;       // O(x,y)

  gal_data_t *data = NULL;
  size_t dsize[2];
  size_t size;
  double *tmp;

  /* Check image type. */
  if (image->type != GAL_TYPE_FLOAT32)
    error (EXIT_FAILURE, 0, "%s: input data must be float 32", __func__);

  /* Process the image and kernel to have same size and be complex numbers. */
  gal_complex_create_padding (image, PSF, &imagepadding, &psfpadding, dsize,
                              &dsize[1]);
  size = dsize[0] * dsize[1]; // Total number of elements.

  /* Normalize and rearange the kernel. */
  gal_complex_normalize (psfpadding, size);
  gal_fft_shift_center (psfpadding, dsize);

  /* Init the solution */
  richardson_lucy_init_solution (&object, size);

  /* Convert to frequency domain. */
  gal_fft_two_dimension_transformation (
      psfpadding, dsize, &psffreq, numthreads, minmapsize, gsl_fft_forward);

  /* Calculate  PSF*(u,v) */
  psffconj = gal_complex_conjugate (psffreq, size);

  for (size_t iteration = 0; iteration < iterations; iteration++)
    {
      richardson_lucy_calculate_next_solution (&object, psffreq, psffconj,
                                               imagepadding, alpha, dsize,
                                               minmapsize, numthreads);
    }

  /* Convert to Real number and convert it to GAL TYPE.*/
  tmp = gal_complex_to_real (object, dsize[0] * dsize[1],
                             COMPLEX_TO_REAL_REAL);
  data = gal_data_alloc (tmp, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1, minmapsize,
                         1, NULL, NULL, NULL); // data has to be 32
  *output = data;

  free (imagepadding);
  free (psfpadding);
  free (psffreq);
  free (psffconj);
  free (object);
}