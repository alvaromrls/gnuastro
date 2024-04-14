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
    // {"multivalue", UI_KEY_MULTIVALUE, "STR", 0,
    //  "This option can take multiple values.", GAL_OPTIONS_GROUP_INPUT,
    //  &p->multivalue, GAL_TYPE_STRLL, GAL_OPTIONS_RANGE_ANY,
    //  GAL_OPTIONS_NOT_MANDATORY, GAL_OPTIONS_NOT_SET,
    //  gal_options_parse_csv_strings_append},

    // {"onoff", UI_KEY_ONOFF, 0, 0, "This option takes no value,.",
    //  GAL_OPTIONS_GROUP_OPERATING_MODE, &p->onoff, GAL_OPTIONS_NO_ARG_TYPE,
    //  GAL_OPTIONS_RANGE_0_OR_1, GAL_OPTIONS_NOT_MANDATORY,
    //  GAL_OPTIONS_NOT_SET},

    {0}};

#endif
