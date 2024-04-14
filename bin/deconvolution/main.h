/*********************************************************************
Deconvolution - A minimal set of files and functions to define a program.
Deconvolution is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef MAIN_H
#define MAIN_H

/* Include necessary headers */
#include <gnuastro/data.h>

#include <gnuastro-internal/options.h>

/* Progarm names.  */
#define PROGRAM_NAME "deconvolution"    /* Program full name.       */
#define PROGRAM_EXEC "astdeconvolution" /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME " (" PACKAGE_NAME ") " PACKAGE_VERSION

/* Main program parameters structure */
struct deconvolution_params {
  /* From command-line */
  struct gal_options_common_params cp; /* Common parameters.           */
  char *filename;   /* Input filename.                         */
  char *kernelname; /* File name of kernel.                    */
  char *khdu;       /* HDU of kernel.                          */
  float lambda;     /*Lambda value.                             */
  /* Output: */
  time_t rawtime; /* Starting time of the program.           */
};

#endif
