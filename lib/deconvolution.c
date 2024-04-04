#include <config.h>

#include <errno.h>
#include <error.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <string.h>

#include <gnuastro/convolve.h>
#include <gnuastro/fits.h>
#include <gnuastro/pointer.h>
#include <gnuastro/threads.h>
#include <gnuastro/tile.h>
#include <gnuastro/wcs.h>

#include <gnuastro-internal/timing.h>

#include <gnuastro/deconvolution.h>

void gal_fft_init(fftparams **params, size_t numthreads, size_t *dim,
                  gsl_const_complex_packed_array input,
                  gsl_complex_packed_array output, gsl_fft_direction sign);
void gal_fft_free(fftparams *params, size_t numthreads);
void gal_fft_one_direction(void *p);

void gal_complex_array_conjugate(gsl_const_complex_packed_array input,
                                 size_t size,
                                 gsl_complex_packed_array *output) {
  gsl_complex_packed_array out;
  out = gal_pointer_allocate(GAL_TYPE_COMPLEX64, size, 1, __func__, "product");
  for (size_t index = 0; index < size; index++) {
    out[index * 2] = input[index * 2];
    out[index * 2 + 1] = -input[index * 2 + 1];
  }
  *output = out;
}

void gal_complex_array_add_scalar(gsl_const_complex_packed_array input,
                                  size_t size, double scalar,
                                  gsl_complex_packed_array *output) {
  gsl_complex_packed_array out;
  out = gal_pointer_allocate(GAL_TYPE_COMPLEX64, size, 1, __func__, "product");
  for (size_t index = 0; index < size; index++) {
    out[index * 2] = input[index * 2] + scalar;
    out[index * 2 + 1] = input[index * 2 + 1];
  }
  *output = out;
}

void gal_fft_complex_array_multiply(gsl_complex_packed_array first,
                                    gsl_complex_packed_array second,
                                    gsl_complex_packed_array *output,
                                    size_t size) {
  gsl_complex_packed_array out;
  out = gal_pointer_allocate(GAL_TYPE_COMPLEX64, size, 1, __func__, "product");

  for (size_t index = 0; index < size; index++) {
    gsl_complex x = {{first[index * 2], first[index * 2 + 1]}};
    gsl_complex y = {{second[index * 2], second[index * 2 + 1]}};
    gsl_complex result = gsl_complex_mul(x, y);
    out[index * 2] = GSL_REAL(result);
    out[index * 2 + 1] = GSL_IMAG(result);
  }
  *output = out;
}

void gal_fft_complex_array_divide(gsl_complex_packed_array first,
                                  gsl_complex_packed_array second,
                                  gsl_complex_packed_array *output,
                                  size_t size) {
  gsl_complex_packed_array out;
  out = gal_pointer_allocate(GAL_TYPE_COMPLEX64, size, 1, __func__, "product");

  for (size_t index = 0; index < size; index++) {
    gsl_complex x = {{first[index * 2], first[index * 2 + 1]}};
    gsl_complex y = {{second[index * 2], second[index * 2 + 1]}};
    if (gsl_complex_abs(y) > MIN_SHARP_SPEC) {
      gsl_complex result = gsl_complex_div(x, y);
      out[index * 2] = GSL_REAL(result);
      out[index * 2 + 1] = GSL_IMAG(result);
    } else {
      out[index * 2] = 0;
      out[index * 2 + 1] = 0;
    }
  }
  *output = out;
}

void helloGNU() { printf("HII ssssssssssssss\n"); }

void gal_create_padding_complex(const gal_data_t *image,
                                const gal_data_t *kernel,
                                padding_complex *output) {
  size_t padsizex, padsizey, totalsize;
  size_t imgsizeX = image->dsize[0];
  size_t imgsizeY = image->dsize[1];
  size_t kernelsizeX = kernel->dsize[0];
  size_t kernelsizeY = kernel->dsize[1];
  float *imagepointer = image->array;
  float *kernelpointer = kernel->array;
  gsl_complex_packed_array pimg;    // pointer to the new image
  gsl_complex_packed_array pkernel; // pointer to the new kernel

  // It's better to have a odd dimension x odd dimension space (so kernel will
  // be centered)
  if (imgsizeX % 2 == 0) {
    padsizex = imgsizeX - 1;
  } else {
    padsizex = imgsizeX;
  }

  if (imgsizeY % 2 == 0) {
    padsizey = imgsizeY - 1;
  } else {
    padsizey = imgsizeY;
  }
  totalsize = padsizex * padsizey;
  // Allocate the image and fill it
  pimg =
      gal_pointer_allocate(GAL_TYPE_COMPLEX64, totalsize, 1, __func__, "pimg");
  for (size_t x = 0; x < padsizex; x++) {
    for (size_t y = 0; y < padsizey; y++) {
      size_t index = imgsizeY * x + y;
      size_t indexcomplex = (padsizey * x + y) * 2;
      pimg[indexcomplex] = (double)imagepointer[index];
      pimg[indexcomplex + 1] = 0.0;
    }
  }
  // Allocate the kernel and fill it
  pkernel = gal_pointer_allocate(GAL_TYPE_COMPLEX64, totalsize, 1, __func__,
                                 "pkernel");
  for (size_t x = 0; x < kernelsizeX; x++) {
    for (size_t y = 0; y < kernelsizeY; y++) {
      size_t kindex = y + x * kernelsizeX;
      size_t index = (x + padsizex / 2 - kernelsizeX / 2) * padsizex +
                     (y + padsizey / 2 - kernelsizeY / 2);
      pkernel[index * 2] = kernelpointer[kindex];
    }
  }
  // Assign values
  output->sizex = padsizex;
  output->sizey = padsizey;
  output->imagepadding = pimg;
  output->kernelpadding = pkernel;
}

void gal_complex_to_real(gsl_complex_packed_array complexarray, size_t size,
                         complex_to_real action, double **output) {
  double *out; // Easier var to access than output
  /* Allocate the space for the real array. */
  out = gal_pointer_allocate(GAL_TYPE_FLOAT64, size, 1, __func__, "out");
  *output = out;

  /* Fill the real array with the derived value from the complex array. */
  switch (action) {
  case COMPLEX_TO_REAL_SPEC:
    for (size_t index = 0; index < size; index++) {
      size_t complexIndex = 2 * index;
      double real = complexarray[complexIndex];
      double im = complexarray[complexIndex + 1];
      out[index] = sqrt((real * real) + (im * im));
    }
    break;
  case COMPLEX_TO_REAL_PHASE:
    for (size_t index = 0; index < size; index++) {
      size_t complexIndex = 2 * index;
      double real = complexarray[complexIndex];
      double im = complexarray[complexIndex + 1];
      out[index] = atan2(real, im);
    }
    break;
  case COMPLEX_TO_REAL_REAL:
    for (size_t index = 0; index < size; index++) {
      size_t complexIndex = 2 * index;
      out[index] = complexarray[complexIndex];
    }
    break;
  default:
    error(EXIT_FAILURE, 0,
          "%s: a bug! Please contact us at %s so we can "
          "correct it. The 'action' code %d is not recognized",
          __func__, PACKAGE_BUGREPORT, action);
  }
}

void gal_normalize_kernel(gsl_complex_packed_array kernel, size_t *dim) {
  gsl_complex_packed_array buffer;
  size_t quadrantxsize = dim[0] / 2;
  size_t quadrantysize = dim[1] / 2;
  // Value normalization
  size_t size = dim[0] * dim[1];
  double cumulativesume = 0.0;
  for (size_t index = 0; index < size; index++) {
    cumulativesume += kernel[index * 2];
  }
  if (cumulativesume == 0.0) {
    error(EXIT_FAILURE, 0,
          "%s: a bug! Please contact us at %s so we can "
          "correct it. Kernel module shouldn't be 0",
          __func__, PACKAGE_BUGREPORT);
  }
  for (size_t index = 0; index < size; index++) {
    kernel[index * 2] /= cumulativesume;
  }
  // Swap diagonal cuadrants for FFT
  buffer =
      gal_pointer_allocate(GAL_TYPE_COMPLEX64, size / 4, 1, __func__, "buffer");
  // 1st and 3rd
  for (size_t x = 0; x < quadrantxsize; x++) {
    for (size_t y = 0; y < quadrantysize; y++) {
      size_t indexb = (x * quadrantxsize + y) * 2;
      size_t index1 = (x * dim[0] + y) * 2;
      size_t index3 = ((x + quadrantxsize) * dim[0] + y + quadrantysize) * 2;
      buffer[indexb] = kernel[index1];
      kernel[index1] = kernel[index3];
      kernel[index3] = buffer[indexb];
    }
  }
  // 2nd and 4th
  for (size_t x = 0; x < quadrantxsize; x++) {
    for (size_t y = 0; y < quadrantysize; y++) {
      size_t indexb = (x * quadrantxsize + y) * 2;
      size_t index2 = (x * dim[0] + y + quadrantysize) * 2;
      size_t index4 = ((x + quadrantxsize) * dim[0] + y) * 2;
      buffer[indexb] = kernel[index2];
      kernel[index2] = kernel[index4];
      kernel[index4] = buffer[indexb];
    }
  }
  free(buffer);
}

void gal_fft_init(fftparams **params, size_t numthreads, size_t *dim,
                  gsl_const_complex_packed_array input,
                  gsl_complex_packed_array output, gsl_fft_direction sign) {
  fftparams *buffer;
  gsl_fft_complex_wavetable *xwave;
  gsl_fft_complex_wavetable *ywave;

  /* Allocate the fftparams array. */
  errno = 0;
  buffer = malloc(numthreads * sizeof(fftparams));
  if (buffer == NULL) {
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for fp", __func__,
          numthreads * sizeof(fftparams));
  }

  xwave = gsl_fft_complex_wavetable_alloc(dim[0]);
  ywave = gsl_fft_complex_wavetable_alloc(dim[1]);

  for (size_t thread = 0; thread < numthreads; thread++) {
    buffer[thread].xwave = xwave;
    buffer[thread].ywave = ywave;
    buffer[thread].xwork = gsl_fft_complex_workspace_alloc(dim[0]);
    buffer[thread].ywork = gsl_fft_complex_workspace_alloc(dim[1]);
    buffer[thread].dim = dim;
    buffer[thread].input = input;
    buffer[thread].output = output;
    buffer[thread].sign = sign;
  }
  *params = buffer;
}

void gal_fft_free(fftparams *params, size_t numthreads) {
  gsl_fft_complex_wavetable_free(params[0].xwave);
  gsl_fft_complex_wavetable_free(params[0].ywave);
  for (size_t thread = 0; thread < numthreads; thread++) {
    gsl_fft_complex_workspace_free(params[thread].xwork);
    gsl_fft_complex_workspace_free(params[thread].ywork);
  }
  free(params);
}

void gal_two_dimension_fft(gsl_const_complex_packed_array input, size_t *dim,
                           gsl_complex_packed_array *output, size_t numthreads,
                           gsl_fft_direction sign) {
  double *out; // Easier var to access than output
  fftparams *params;
  size_t size = dim[0] * dim[1];
  char *mmapname = NULL;
  size_t *thrds, thrdscols;
  /* Allocate the space for the real array. */
  printf("Allocating fourier img\n");
  out = gal_pointer_allocate(GAL_TYPE_COMPLEX64, size, 1, __func__, "outf");
  printf("init fft params\n");
  memcpy(out, input, size * sizeof(gsl_complex_packed_array) * 2);
  *output = out;
  gal_fft_init(&params, numthreads, dim, input, out, sign);

  /* 1D FFT on each row. */
  mmapname = gal_threads_dist_in_threads(dim[0], numthreads, MIN_MAP_SIZE, 0,
                                         &thrds, &thrdscols);
  if (numthreads == 1) {
    params[0].stride = 1;
    params[0].indexs = thrds;
    gal_fft_one_direction(params);
  } else {
    // todo
  }

  /* Clean up. */
  if (mmapname) {
    gal_pointer_mmap_free(&mmapname, 0);
  } else {
    free(thrds);
  }

  /* 1D FFT on each column. */
  mmapname = gal_threads_dist_in_threads(dim[1], numthreads, MIN_MAP_SIZE, 0,
                                         &thrds, &thrdscols);
  if (numthreads == 1) {
    params[0].stride = dim[0];
    params[0].indexs = thrds;
    gal_fft_one_direction(params);
  } else {
    // todo
  }
  /* Clean up. */
  if (mmapname) {
    gal_pointer_mmap_free(&mmapname, 0);
  } else {
    free(thrds);
  }
  printf("clean fft params\n");
  gal_fft_free(params, numthreads);
}

void gal_fft_one_direction(void *p) {
  fftparams *params = (fftparams *)p;
  size_t lenght;
  size_t indmultip;
  size_t *indexs = params->indexs;
  gsl_fft_complex_workspace *work;
  gsl_fft_complex_wavetable *wavetable;
  gsl_complex_packed_array data;
  if (params->stride == 1) {
    lenght = params->dim[1];
    work = params->xwork;
    wavetable = params->xwave;
    indmultip = params->dim[1];
  } else {
    lenght = params->dim[0];
    work = params->xwork;
    wavetable = params->xwave;
    indmultip = 1;
  }
  for (size_t i = 0; indexs[i] != GAL_BLANK_SIZE_T; ++i) {
    data = &params->output[2 * indexs[i] * indmultip];
    gsl_fft_complex_transform(data, params->stride, lenght, wavetable, work,
                              params->sign);
  }
}

void gal_deconvolution_tikhonov(const gal_data_t *image, const gal_data_t *PSF,
                                double lambda, gal_data_t **output) {
  padding_complex d;
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
  gal_create_padding_complex(image, PSF, &d);
  dsize[0] = d.sizex;
  dsize[1] = d.sizey;
  size = dsize[0] * dsize[1];
  gal_normalize_kernel(d.kernelpadding, dsize);
  gal_two_dimension_fft(d.kernelpadding, dsize, &psffreq, 1, gsl_fft_forward);

  gal_two_dimension_fft(d.imagepadding, dsize, &imagefreq, 1, gsl_fft_forward);
  gal_complex_array_conjugate(psffreq, size, &psffconj);
  gal_fft_complex_array_multiply(psffreq, psffconj, &psffreqsquare, size);
  gal_fft_complex_array_multiply(psffreq, imagefreq, &numerator, size);
  gal_complex_array_add_scalar(psffreqsquare, size, lambda, &denominator);
  gal_fft_complex_array_divide(numerator, denominator, &deconvolutionfreq,
                               size);
  gal_two_dimension_fft(deconvolutionfreq, dsize, &deconvolution, 1,
                        gsl_fft_backward);
  gal_complex_to_real(deconvolution, d.sizex * d.sizey, COMPLEX_TO_REAL_REAL,
                      &tmp);
  data = gal_data_alloc(tmp, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1, MIN_MAP_SIZE,
                        1, NULL, NULL, NULL);
  *output = data;
  free(d.imagepadding);
  free(d.kernelpadding);
  free(psffreq);
  free(psffconj);
  free(imagefreq);
  free(psffreqsquare);
  free(numerator);
  free(denominator);
  free(deconvolutionfreq);
  free(deconvolution);
}
