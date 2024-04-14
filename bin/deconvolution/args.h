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
#ifndef ARGS_H
#define ARGS_H

/* Array of acceptable options. */
struct argp_option program_options[] = {
  /* Inputs */
  { "kernel", UI_KEY_KERNEL, "FITS", 0, "File name of kernel for convolution.",
    GAL_OPTIONS_GROUP_INPUT, &p->kernelname, GAL_TYPE_STRING,
    GAL_OPTIONS_RANGE_ANY, GAL_OPTIONS_NOT_MANDATORY, GAL_OPTIONS_NOT_SET },
  { "khdu", UI_KEY_KHDU, "STR", 0, "HDU containing the kernel.",
    GAL_OPTIONS_GROUP_INPUT, &p->khdu, GAL_TYPE_STRING, GAL_OPTIONS_RANGE_ANY,
    GAL_OPTIONS_MANDATORY, GAL_OPTIONS_NOT_SET },
  { "lambda", UI_KEY_LAMBDA, "FLT", 0, "Deconvolution: value for Tikhonov alg",
    GAL_OPTIONS_GROUP_INPUT, &p->lambda, GAL_TYPE_FLOAT64,
    GAL_OPTIONS_RANGE_GE_0_LE_1, GAL_OPTIONS_MANDATORY, GAL_OPTIONS_NOT_SET },
  { "algorithm", UI_KEY_ALGORITHM, "STR", 0,
    "Deconvolution: Algorithm selection", GAL_OPTIONS_GROUP_INPUT,
    &p->algorithmstr, GAL_TYPE_STRING, GAL_OPTIONS_RANGE_ANY,
    GAL_OPTIONS_MANDATORY, GAL_OPTIONS_NOT_SET },
  { 0 }
};

/* Define the child argp structure. */
struct argp gal_options_common_child = { gal_commonopts_options,
                                         gal_options_common_argp_parse,
                                         NULL,
                                         NULL,
                                         NULL,
                                         NULL,
                                         NULL };

/* Use the child argp structure in list of children (only one for now). */
struct argp_child children[]
    = { { &gal_options_common_child, 0, NULL, 0 }, { 0, 0, 0, 0 } };

/* Set all the necessary argp parameters. */
struct argp thisargp
    = { program_options, parse_opt, args_doc, doc, children, NULL, NULL };
#endif
