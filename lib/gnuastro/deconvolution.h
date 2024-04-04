/*********************************************************************
Functions to interface with DS9 files.
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
#ifndef __GAL_DECONVOLUTION_H__
#define __GAL_DECONVOLUTION_H__

#include <gnuastro/data.h>

#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN_MAP_SIZE 9223372036854775807UL
#define MIN_SHARP_SPEC 0.0000000050
typedef struct {
  size_t sizex;
  size_t sizey;
  gsl_complex_packed_array imagepadding;
  gsl_complex_packed_array kernelpadding;
} padding_complex;

typedef enum {
  COMPLEX_TO_REAL_INVALID, /* ==0 by C standard. */

  COMPLEX_TO_REAL_SPEC,
  COMPLEX_TO_REAL_PHASE,
  COMPLEX_TO_REAL_REAL,
} complex_to_real;

void helloGNU();

/**
 * @brief Function to generate a padding in kernel (and image if needed) and
 * converts them from real images to complex ones
 *
 * @param image
 * @param kernel
 * @param output
 */
void gal_create_padding_complex(const gal_data_t *image,
                                const gal_data_t *kernel,
                                padding_complex *output);

void gal_complex_to_real(gsl_complex_packed_array complexarray, size_t size,
                         complex_to_real action, double **output);

/**
 * @brief Normalize the real part from a given kernel and swaps quadrants 1-3
 * and 2-4. Needed for fft transformations.
 * 2 | 3
 * ------
 * 1 | 4
 *
 * @param kernel
 * @param size
 */
void gal_normalize_kernel(gsl_complex_packed_array kernel, size_t *dim);

typedef struct {
  size_t id;
  size_t stride;
  gsl_fft_complex_wavetable *xwave;
  gsl_fft_complex_wavetable *ywave;
  gsl_fft_complex_workspace *xwork;
  gsl_fft_complex_workspace *ywork;
  size_t *indexs;
  pthread_barrier_t *b;
  gsl_const_complex_packed_array input;
  size_t *dim;
  gsl_complex_packed_array output;
  gsl_fft_direction sign;
} fftparams;

void gal_two_dimension_fft(gsl_const_complex_packed_array input, size_t *dim,
                           gsl_complex_packed_array *output, size_t numthreads,
                           gsl_fft_direction sign);

void gal_fft_complex_array_multiply(gsl_complex_packed_array first,
                                    gsl_complex_packed_array second,
                                    gsl_complex_packed_array *output,
                                    size_t size);

void gal_fft_complex_array_divide(gsl_complex_packed_array first,
                                  gsl_complex_packed_array second,
                                  gsl_complex_packed_array *output,
                                  size_t size);

void gal_complex_array_conjugate(gsl_const_complex_packed_array input,
                                 size_t size, gsl_complex_packed_array *output);

void gal_complex_array_add_scalar(gsl_const_complex_packed_array input,
                                  size_t size, double scalar,
                                  gsl_complex_packed_array *output);

void gal_deconvolution_tikhonov(const gal_data_t *image, const gal_data_t *PSF,
                                double lambda, gal_data_t **output);

#endif