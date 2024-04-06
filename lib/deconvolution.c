#include <config.h>

#include <errno.h>
#include <error.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <string.h>

#include <gnuastro/complexarray.h>
#include <gnuastro/convolve.h>
#include <gnuastro/fft.h>
#include <gnuastro/fits.h>
#include <gnuastro/pointer.h>
#include <gnuastro/tile.h>
#include <gnuastro/wcs.h>

#include <gnuastro-internal/timing.h>

#include <gnuastro/deconvolution.h>

void
gal_deconvolution_tikhonov (const gal_data_t *image, const gal_data_t *PSF,
                            double lambda, gal_data_t **output)
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
  gsl_complex_packed_array deconvolution;
  gal_data_t *data = NULL;
  size_t dsize[2];
  size_t size;
  double *tmp;

  if (image->type != GAL_TYPE_FLOAT32)
    error (EXIT_FAILURE, 0, "%s: input data must be float 64", __func__);

  gal_complex_array_create_padding (image, PSF, &imagepadding, &psfpadding,
                                    dsize, &dsize[1]);
  size = dsize[0] * dsize[1];
  gal_complex_array_normalize (psfpadding, size);
  gal_fft_swap_quadrant (psfpadding, dsize);
  gal_fft_two_dimension_transformation (psfpadding, dsize, &psffreq, 1,
                                        gsl_fft_forward);
  gal_fft_two_dimension_transformation (imagepadding, dsize, &imagefreq, 1,
                                        gsl_fft_forward);
  gal_complex_array_conjugate (psffreq, size, &psffconj);
  gal_complex_array_multiply (psffreq, psffconj, &psffreqsquare, size);
  gal_complex_array_multiply (psffreq, imagefreq, &numerator, size);
  gal_complex_array_add_scalar (psffreqsquare, size, lambda + I * 0,
                                &denominator);
  gal_complex_array_divide (numerator, denominator, &deconvolutionfreq, size,
                            lambda);
  gal_fft_two_dimension_transformation (deconvolutionfreq, dsize,
                                        &deconvolution, 1, gsl_fft_backward);
  gal_complex_array_to_real (deconvolution, dsize[0] * dsize[1],
                             COMPLEX_TO_REAL_REAL, &tmp);
  data
      = gal_data_alloc (tmp, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1, MIN_MAP_SIZE,
                        1, NULL, NULL, NULL); // data has to be 32
  *output = data;
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
