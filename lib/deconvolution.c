#include <config.h>

#include <errno.h>
#include <error.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <string.h>

#include <gnuastro/complexarray.h>
#include <gnuastro/convolve.h>
#include <gnuastro/fits.h>
#include <gnuastro/pointer.h>
#include <gnuastro/threads.h>
#include <gnuastro/tile.h>
#include <gnuastro/wcs.h>

#include <gnuastro-internal/timing.h>

#include <gnuastro/deconvolution.h>

void gal_fft_init (fftparams **params, size_t numthreads, size_t *dim,
                   gsl_const_complex_packed_array input,
                   gsl_complex_packed_array output, gsl_fft_direction sign);
void gal_fft_free (fftparams *params, size_t numthreads);
void gal_fft_one_direction (void *p);

void
gal_swap_quadrant (gsl_complex_packed_array kernel, size_t *dim)
{
  gsl_complex_packed_array buffer;
  size_t quadrantxsize = dim[0] / 2;
  size_t quadrantysize = dim[1] / 2;
  size_t size = dim[0] * dim[1];
  // Swap diagonal cuadrants for FFT
  buffer = gal_pointer_allocate (GAL_TYPE_COMPLEX64, size / 4, 1, __func__,
                                 "buffer");
  // 1st and 3rd
  for (size_t x = 0; x < quadrantxsize; x++)
    {
      for (size_t y = 0; y < quadrantysize; y++)
        {
          size_t indexb = (x * quadrantxsize + y) * 2;
          size_t index1 = (x * dim[0] + y) * 2;
          size_t index3
              = ((x + quadrantxsize) * dim[0] + y + quadrantysize) * 2;
          buffer[indexb] = kernel[index1];
          kernel[index1] = kernel[index3];
          kernel[index3] = buffer[indexb];
        }
    }
  // 2nd and 4th
  for (size_t x = 0; x < quadrantxsize; x++)
    {
      for (size_t y = 0; y < quadrantysize; y++)
        {
          size_t indexb = (x * quadrantxsize + y) * 2;
          size_t index2 = (x * dim[0] + y + quadrantysize) * 2;
          size_t index4 = ((x + quadrantxsize) * dim[0] + y) * 2;
          buffer[indexb] = kernel[index2];
          kernel[index2] = kernel[index4];
          kernel[index4] = buffer[indexb];
        }
    }
  free (buffer);
}

void
gal_fft_init (fftparams **params, size_t numthreads, size_t *dim,
              gsl_const_complex_packed_array input,
              gsl_complex_packed_array output, gsl_fft_direction sign)
{
  fftparams *buffer;
  gsl_fft_complex_wavetable *xwave;
  gsl_fft_complex_wavetable *ywave;

  /* Allocate the fftparams array. */
  errno = 0;
  buffer = malloc (numthreads * sizeof (fftparams));
  if (buffer == NULL)
    {
      error (EXIT_FAILURE, errno, "%s: allocating %zu bytes for fp", __func__,
             numthreads * sizeof (fftparams));
    }

  xwave = gsl_fft_complex_wavetable_alloc (dim[0]);
  ywave = gsl_fft_complex_wavetable_alloc (dim[1]);

  for (size_t thread = 0; thread < numthreads; thread++)
    {
      buffer[thread].xwave = xwave;
      buffer[thread].ywave = ywave;
      buffer[thread].xwork = gsl_fft_complex_workspace_alloc (dim[0]);
      buffer[thread].ywork = gsl_fft_complex_workspace_alloc (dim[1]);
      buffer[thread].dim = dim;
      buffer[thread].input = input;
      buffer[thread].output = output;
      buffer[thread].sign = sign;
    }
  *params = buffer;
}

void
gal_fft_free (fftparams *params, size_t numthreads)
{
  gsl_fft_complex_wavetable_free (params[0].xwave);
  gsl_fft_complex_wavetable_free (params[0].ywave);
  for (size_t thread = 0; thread < numthreads; thread++)
    {
      gsl_fft_complex_workspace_free (params[thread].xwork);
      gsl_fft_complex_workspace_free (params[thread].ywork);
    }
  free (params);
}

void
gal_two_dimension_fft (gsl_const_complex_packed_array input, size_t *dim,
                       gsl_complex_packed_array *output, size_t numthreads,
                       gsl_fft_direction sign)
{
  double *out; // Easier var to access than output
  fftparams *params;
  size_t size = dim[0] * dim[1];
  char *mmapname = NULL;
  size_t *thrds, thrdscols;
  /* Allocate the space for the real array. */
  printf ("Allocating fourier img\n");
  out = gal_pointer_allocate (GAL_TYPE_COMPLEX64, size, 1, __func__, "outf");
  printf ("init fft params\n");
  memcpy (out, input, size * sizeof (gsl_complex_packed_array) * 2);
  gal_fft_init (&params, numthreads, dim, input, out, sign);

  /* 1D FFT on each row. */
  mmapname = gal_threads_dist_in_threads (dim[0], numthreads, MIN_MAP_SIZE, 0,
                                          &thrds, &thrdscols);
  if (numthreads == 1)
    {
      params[0].stride = 1;
      params[0].indexs = thrds;
      gal_fft_one_direction (params);
    }
  else
    {
      // todo
    }

  /* Clean up. */
  if (mmapname)
    {
      gal_pointer_mmap_free (&mmapname, 0);
    }
  else
    {
      free (thrds);
    }

  /* 1D FFT on each column. */
  mmapname = gal_threads_dist_in_threads (dim[1], numthreads, MIN_MAP_SIZE, 0,
                                          &thrds, &thrdscols);
  if (numthreads == 1)
    {
      params[0].stride = dim[0];
      params[0].indexs = thrds;
      gal_fft_one_direction (params);
    }
  else
    {
      // todo
    }
  /* Clean up. */
  if (mmapname)
    {
      gal_pointer_mmap_free (&mmapname, 0);
    }
  else
    {
      free (thrds);
    }
  if (sign == gsl_fft_backward)
    {
      gal_complex_array_scale (out, (1.0 / (double)size), size);
    }
  printf ("clean fft params\n");
  gal_fft_free (params, numthreads);
  *output = out;
}

void
gal_fft_one_direction (void *p)
{
  fftparams *params = (fftparams *)p;
  size_t lenght;
  size_t indmultip;
  size_t *indexs = params->indexs;
  gsl_fft_complex_workspace *work;
  gsl_fft_complex_wavetable *wavetable;
  gsl_complex_packed_array data;
  if (params->stride == 1)
    {
      lenght = params->dim[1];
      work = params->xwork;
      wavetable = params->xwave;
      indmultip = params->dim[1];
    }
  else
    {
      lenght = params->dim[0];
      work = params->xwork;
      wavetable = params->xwave;
      indmultip = 1;
    }
  for (size_t i = 0; indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      data = &params->output[2 * indexs[i] * indmultip];
      gsl_fft_complex_transform (data, params->stride, lenght, wavetable, work,
                                 params->sign);
    }
}

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
  gal_swap_quadrant (psfpadding, dsize);
  gal_two_dimension_fft (psfpadding, dsize, &psffreq, 1, gsl_fft_forward);
  gal_two_dimension_fft (imagepadding, dsize, &imagefreq, 1, gsl_fft_forward);
  gal_complex_array_conjugate (psffreq, size, &psffconj);
  gal_complex_array_multiply (psffreq, psffconj, &psffreqsquare, size);
  gal_complex_array_multiply (psffreq, imagefreq, &numerator, size);
  gal_complex_array_add_scalar (psffreqsquare, size, lambda, &denominator);
  gal_complex_array_divide (numerator, denominator, &deconvolutionfreq, size,
                            lambda + I * 0);
  gal_two_dimension_fft (deconvolutionfreq, dsize, &deconvolution, 1,
                         gsl_fft_backward);
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
