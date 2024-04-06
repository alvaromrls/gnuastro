/*********************************************************************
Convolve - Convolve input data with a given kernel.
Convolve is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2024 Free Software Foundation, Inc.

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

#include <stdio.h>
#include <stdlib.h>

#include <gnuastro/deconvolution.h>

#include <gnuastro-internal/timing.h>

#include "main.h"

#include "convolve.h"
#include "ui.h" /* Needs convolveparams in main.h */

int
main (int argc, char *argv[])
{
  struct timeval t1;
  struct convolveparams p = { { { 0 }, 0 }, 0 };
  gal_data_t *data = NULL;
  /* Set the starting time. */
  time (&p.rawtime);
  gettimeofday (&t1, NULL);
  /* Read the input parameters. */
  ui_read_check_inputs_setup (argc, argv, &p);
  /* Run Convolve. */
  // convolve(&p);
  /* Save the padded kernel image. */
  /* Prepare the data structure for viewing the steps, note that we
    don't need the array that is initially made. */
  // data = gal_data_alloc (NULL, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1,
  //                        p.cp.minmapsize, p.cp.quietmmap, NULL, NULL, NULL);
  // gal_complex_to_real (d.kernelpadding, d.sizex * d.sizey,
  //                      COMPLEX_TO_REAL_REAL, &tmp);
  // data->array = tmp;
  // gal_fits_img_write (data, "kernelpadeado.fits", NULL, 0);
  // gal_complex_to_real (d.kernelpadding, d.sizex * d.sizey,
  //                      COMPLEX_TO_REAL_REAL, &tmp);
  // data = gal_data_alloc (tmp, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1,
  gal_deconvolution_tikhonov (p.input, p.kernel, 100, p.cp.numthreads, &data);

  gal_fits_img_write (data, "deconvolution.fits", NULL, 0);
  free (data);

  /* Free all non-freed allocations. */
  ui_free_report (&p, &t1);

  /* Return successfully. */
  return EXIT_SUCCESS;
}
