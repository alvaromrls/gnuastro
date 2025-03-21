/*********************************************************************
Match - A program to match catalogs and WCS warps
Match is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2025 Free Software Foundation, Inc.

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

#include <argp.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <string.h>

#include <gnuastro/fits.h>
#include <gnuastro/match.h>
#include <gnuastro/threads.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/options.h>
#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/fixedstringmacros.h>

#include "main.h"

#include "ui.h"
#include "authors-cite.h"





/**************************************************************/
/*********      Argp necessary global entities     ************/
/**************************************************************/
/* Definition parameters for the Argp: */
const char *
argp_program_version = PROGRAM_STRING "\n"
                       GAL_STRINGS_COPYRIGHT
                       "\n\nWritten/developed by "PROGRAM_AUTHORS;

const char *
argp_program_bug_address = PACKAGE_BUGREPORT;

static char
args_doc[] = "ASTRdata";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" matches catalogs of objects "
  "and (by default) will return the re-arranged matching inputs. The "
  "optional log file will return low-level information about the match "
  "(indexs and distances).\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct matchparams *p,
                      struct argp_option *program_options,
                      struct argp_option *gal_commonopts_options)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;


  /* Set the necessary common parameters structure. */
  cp->program_struct     = p;
  cp->poptions           = program_options;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->numthreads         = gal_threads_number();
  cp->coptions           = gal_commonopts_options;


  /* Modify common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    {
      /* Select individually. */
      switch(cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_HDU:
          cp->coptions[i].doc="Extension name or number of first input.";
          break;
        case GAL_OPTIONS_KEY_TYPE:
          cp->coptions[i].flags=OPTION_HIDDEN;
          break;
        }

      /* Select by group. */
      switch(cp->coptions[i].group)
        {
        case GAL_OPTIONS_GROUP_TESSELLATION:
          cp->coptions[i].doc=NULL; /* Necessary to remove title. */
          cp->coptions[i].flags=OPTION_HIDDEN;
          break;
        }
    }
}





/* Parse a single option: */
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct matchparams *p = state->input;

  /* Pass 'gal_options_common_params' into the child parser.  */
  state->child_inputs[0] = &p->cp;

  /* In case the user incorrectly uses the equal sign (for example
     with a short format or with space in the long format, then 'arg'
     start with (if the short version was called) or be (if the long
     version was called with a space) the equal sign. So, here we
     check if the first character of arg is the equal sign, then the
     user is warned and the program is stopped: */
  if(arg && arg[0]=='=')
    argp_error(state, "incorrect use of the equal sign ('='). For short "
               "options, '=' should not be used and for long options, "
               "there should be no space between the option, equal sign "
               "and value");

  /* Set the key to this option. */
  switch(key)
    {
    /* Read the non-option tokens (arguments). */
    case ARGP_KEY_ARG:
      /* The user may give a shell variable that is empty! In that case
         'arg' will be an empty string! We don't want to account for such
         cases (and give a clear error that no input has been given). */
      if(p->input1name)
        {
          if(p->input2name)
            argp_error(state, "only two arguments (input files) "
                       "should be given, not any more");
          else
            { if(arg[0]!='\0') p->input2name=arg; }
        }
      else
        { if(arg[0]!='\0') p->input1name=arg; }
      break;


    /* This is an option, set its value. */
    default:
      return gal_options_set_from_key(key, arg, p->cp.poptions, &p->cp);
    }

  return 0;
}





void *
ui_parse_match_type(struct argp_option *option, char *arg,
                    char *filename, size_t lineno, void *junk)
{
  char *str=NULL;
  uint8_t value=GAL_MATCH_ARRANGE_INVALID;

  /* If we are printing values. */
  if(lineno==-1)
    {
      /* The output must be an allocated string (will be 'free'd later). */
      value=*(uint8_t *)(option->value);
      switch(value)
        {
        case GAL_MATCH_ARRANGE_FULL:
          gal_checkset_allocate_copy("full", &str);
          break;
        case GAL_MATCH_ARRANGE_INNER:
          gal_checkset_allocate_copy("inner", &str);
          break;
        case GAL_MATCH_ARRANGE_OUTER:
          gal_checkset_allocate_copy("outer", &str);
          break;
        case GAL_MATCH_ARRANGE_OUTERWITHINAPERTURE:
          gal_checkset_allocate_copy("outer-within-aperture", &str);
          break;
        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' "
                "to fix the problem. %u is not a recognized match type",
                __func__, PACKAGE_BUGREPORT, value);
        }
      return str;
    }

  /* If we are reading/parsing values. */
  else
    {
      /* If the option is already set, just return. */
      if(option->set) return NULL;

      /* Read the value. */
      if(      !strcmp(arg, "full") )  value = GAL_MATCH_ARRANGE_FULL;
      else if( !strcmp(arg, "inner") ) value = GAL_MATCH_ARRANGE_INNER;
      else if( !strcmp(arg, "outer") ) value = GAL_MATCH_ARRANGE_OUTER;
      else if( !strcmp(arg, "outer-within-aperture") )
        value = GAL_MATCH_ARRANGE_OUTERWITHINAPERTURE;
      else
        error_at_line(EXIT_FAILURE, 0, filename, lineno, "'%s' (value "
                      "to '%s' option) couldn't be recognized as a known "
                      "mathing type. Acceptable values are 'full', "
                      "'inner', 'outer' or 'outer-within-aperture'", arg,
                      option->name);
      *(uint8_t *)(option->value)=value;

      /* For no un-used variable warning. This function doesn't need the
         pointer. */
      return junk=NULL;
    }
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* Check ONLY the options. When arguments are involved, do the check
   in 'ui_check_options_and_arguments'. */
static void
ui_check_only_options(struct matchparams *p)
{
  /* Set the k-d tree mode. */
  if(      !strcmp(p->kdtree,"build")    )
    p->kdtreemode=MATCH_KDTREE_BUILD;
  else if( !strcmp(p->kdtree,"internal") )
    p->kdtreemode=MATCH_KDTREE_INTERNAL;
  else if( !strcmp(p->kdtree,"disable")  )
    p->kdtreemode=MATCH_KDTREE_DISABLE;
  else if( gal_fits_name_is_fits(p->kdtree) )
    p->kdtreemode=MATCH_KDTREE_FILE;
  else
    error(EXIT_FAILURE, 0, "'%s' is not valid for '--kdtree'. The "
          "following values are accepted: 'build' (to build the k-d "
          "tree in the file given to '--output'), 'internal' (to "
          "force internal usage of a k-d tree for the matching), "
          "'disable' (to not use a k-d tree at all), a FITS file "
          "name (the file to read a created k-d tree from)", p->kdtree);

  /* Make sure that the k-d tree build mode is not called with
     '--outcols'. */
  if( p->kdtreemode==MATCH_KDTREE_BUILD && (p->outcols || p->coord) )
    error(EXIT_FAILURE, 0, "the '--kdtree=build' option is "
          "incompatible with the '--outcols' or '--coord' options "
          "(because in the k-d tree building mode doesn't involve "
          "actual matching. It will only build k-d tree and write "
          "it to a file so it can be used in future matches)");

  /* Make sure that a HDU is also specified for the k-d tree when its an
     external file. */
  if( p->kdtreemode==MATCH_KDTREE_FILE && p->kdtreehdu==NULL )
    error(EXIT_FAILURE, 0, "no '--kdtreehdu' specified! When a FITS "
          "file is given to '--kdtree', you should specify its "
          "extension/HDU within the file using '--kdtreehdu'. You "
          "can either use the HDU name, or its number (counting "
          "from 0). If you aren't sure what HDUs exist in the file, "
          "you can use the 'astfits %s' command to see the full list",
          p->kdtree);

  /* For an outer match, we only operate in k-d tree mode. */
  if( (   p->type==GAL_MATCH_ARRANGE_OUTER
       || p->type==GAL_MATCH_ARRANGE_OUTERWITHINAPERTURE )
      && p->kdtreemode==MATCH_KDTREE_DISABLE )
    error(EXIT_FAILURE, 0, "the 'outer' arrangement only works with the "
          "k-d tree algorithm (in other words, it conflicts with "
          "'--kdtree=disable' option). Please set '--kdtree=internal' "
          "(if you want to build the k-d tree very time), or to "
          "'--kdtree=file.fits' (where 'file.fits' is an existing k-d "
          "tree that was built once and can be used many times for "
          "a faster result on a repeated reference catalog). See the "
          "description of the '--kdtree' option of the Match program "
          "in Gnuastro manual for more (by running 'info %s')",
          PROGRAM_EXEC);

  /* The '--notmatched' option is only valid for the inner type of
     matching. */
  if(p->notmatched && p->type!=GAL_MATCH_ARRANGE_INNER)
    error(EXIT_FAILURE, 0, "the '--notmatched' option is only defined "
          "with the '--arrange=inner' match type");
}





/* Two necessary catalogs: First:  Standard input, or a file.
                           Second: '--coord',      or a file.
 */
static void
ui_check_options_and_arguments(struct matchparams *p)
{
  /* When '--coord' is given, or we are in k-d tree building mode, there
     should be no second catalog (argument). */
  if(p->coord || p->kdtreemode==MATCH_KDTREE_BUILD)
    {
      /* Make sure no second argument is given. */
      if(p->input2name)
        error(EXIT_FAILURE, 0, "only one argument can be given with the "
              "'--coord' or '--kdtree=build' options");

      /* In case '--ccol2' is given. */
      if(p->ccol2 && p->cp.quiet==0)
        {
          /* Print a warning to let the user know the option will not be
             used (the user had probably confused things if they have given
             it). */
          error(EXIT_SUCCESS, 0, "WARNING: with '--coord' or '--kdtree', "
                "you only need to set '--ccol1' (since there is only a "
                "single input table), '--ccol2' is not needed");
        }
    }

  /* We need two input catalogs. */
  else
    {
      /* Without 'coord' atleast one input is necessary. */
      if(p->input1name==NULL)
        error(EXIT_FAILURE, 0, "no inputs!\n\n"
              "Two inputs are necessary. The first can be a file, or from "
              "the standard input (e.g., a pipe). The second can be a "
              "file, or its coordinates can be directly specified on the "
              "command-line with '--coord'");

      /* If the first input should be read from the standard input, the
         contents of 'input1name' actually belong to 'input2name'. */
      if(p->input2name==NULL)
        {
          p->input2name=p->input1name;
          p->input1name=NULL;
        }
    }

  /* If the first input is a FITS file, make sure its HDU is also given. */
  if( p->input1name
      && gal_fits_file_recognized(p->input1name)
      && p->cp.hdu==NULL )
    error(EXIT_FAILURE, 0, "no HDU for first input. When the input is "
          "a FITS file, a HDU must also be specified, you can use the "
          "'--hdu' ('-h') option and give it the HDU number (starting "
          "from zero), extension name, or anything acceptable by "
          "CFITSIO");
  if( p->input2name
      && gal_fits_file_recognized(p->input2name)
      && p->hdu2==NULL )
    error(EXIT_FAILURE, 0, "no HDU for second input. Please use "
          "the '--hdu2' ('-H') option and give it the HDU number "
          "(starting from zero), extension name, or anything "
          "acceptable by CFITSIO");
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
/* The final aperture must have the following values:

       p->aperture[0]: Major axis length.
       p->aperture[1]: Axis ratio.
       p->aperture[2]: Position angle (relative to first dim).     */
static void
ui_read_columns_aperture_2d(struct matchparams *p)
{
  size_t apersize=3;
  gal_data_t *newaper=NULL;
  double *naper=NULL, *oaper=p->aperture->array;

  /* A general sanity check: the first two elements of aperture cannot be
     zero or negative. */
  if( oaper[0]<=0 )
    error(EXIT_FAILURE, 0, "the first value of '--aperture' cannot be "
          "zero or negative");
  if( p->aperture->size>1 && oaper[1]<=0 )
    error(EXIT_FAILURE, 0, "the second value of '--aperture' cannot be "
          "zero or negative");

  /* Will be needed in more than one case. */
  if(p->aperture->size!=apersize)
    {
      newaper=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &apersize, NULL,
                             0, -1, 1, NULL, NULL, NULL);
      naper=newaper->array;
    }

  /* Different based on  */
  switch(p->aperture->size)
    {
    case 1:
      naper[0]=oaper[0];
      naper[1]=1;
      naper[2]=0;
      break;

    case 2:
      naper[0] = oaper[0]>oaper[1] ? oaper[0]          : oaper[1];
      naper[1] = oaper[0]>oaper[1] ? oaper[1]/oaper[0] : oaper[0]/oaper[1];
      naper[2] = oaper[0]>oaper[1] ? 0                 : 90;
      break;

    case 3:
      if(oaper[1]>1)
        error(EXIT_FAILURE, 0, "second value to '--aperture' is larger "
              "than one. When three numbers are given to this option, "
              "the second is the axis ratio (which must always be less "
              "than 1).");
      break;

    default:
      error(EXIT_FAILURE, 0, "%zu values given to '--aperture'. In 2D, "
            "this option can only take 1, 2, or 3 values",
            p->aperture->size);
    }

  /* If a new aperture was defined, then replace it with the exitsting
     one. */
  if(newaper)
    {
      gal_data_free(p->aperture);
      p->aperture=newaper;
    }
}





/* The 3D aperture must finally have these elements.

      p->aperture[0] = Major axis length.
      p->aperture[1] = Second axis ratio.
      p->aperture[2] = Third axis ratio.
      p->aperture[3] = First ZXZ Euler angle rotation.
      p->aperture[4] = Second ZXZ Euler angle rotation.
      p->aperture[5] = Third ZXZ Euler angle rotation.   */
static void
ui_read_columns_aperture_3d(struct matchparams *p)
{
  size_t apersize=6;
  gal_data_t *newaper=NULL;
  double *naper=NULL, *oaper=p->aperture->array;

  /* A general sanity check: the first two elements of aperture cannot be
     zero or negative. */
  if( oaper[0]<=0 )
    error(EXIT_FAILURE, 0, "the first value of '--aperture' cannot be "
          "zero or negative");
  if( p->aperture->size>2 && (oaper[1]<=0 || oaper[2]<=0) )
    error(EXIT_FAILURE, 0, "the second and third values of '--aperture' "
          "cannot be zero or negative");

  /* Will be needed in more than one case. */
  if(p->aperture->size!=apersize)
    {
      newaper=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &apersize, NULL,
                             0, -1, 1, NULL, NULL, NULL);
      naper=newaper->array;
    }

  /* Different based on  */
  switch(p->aperture->size)
    {
    case 1:
      naper[0] = oaper[0];
      naper[1] = naper[2] = 1;
      naper[3] = naper[4] = naper[5] = 0;
      break;

    case 3:
      /* Major axis is along the first dimension, no rotation necessary. */
      if(oaper[0]>=oaper[1] && oaper[0]>=oaper[2])
        {
          naper[0] = oaper[0];
          naper[1] = oaper[1] / oaper[0];
          naper[2] = oaper[2] / oaper[0];
          naper[3] = naper[4] = naper[5] = 0;
        }

      /* Major axis is along the second dimension. So we want 'X' to be in
         the direction of 'y'. Therefore, just the first Eurler ZXZ
         rotation is necessary by 90 degrees. Here is how the rotated
         coordinates (X,Y,Z) look like (after one rotation about Z).

                  |z (before)                   Z|  /Y
                  |                              | /
                  |                  ==>         |/
                  .-------- y                    .------X
                 /
              x /

         You see how the major axis (X) now lies in the second original
         direction (y). The length along 'x' is now along 'Y' and the third
         hasn't changed. Note that we are talking about the semi-axis
         lengths (a scalar, not a vector), so +- is irrelevant. */
      else if(oaper[1]>=oaper[0] && oaper[1]>=oaper[2])
        {
          naper[0] = oaper[1];
          naper[1] = oaper[0] / oaper[1];
          naper[2] = oaper[2] / oaper[1];
          naper[3] = 90;
          naper[4] = naper[5] = 0;
        }

      /* The major axis is along the third dimension. So we want 'X' to
         point in the direction of 'z'. To get to that point, we need 90
         degree rotations in all three Euler ZXZ rotations.

              |z (before)     z'|  /y'            |y''                  |X
              |                 | /               |                     |
              |            ==>  |/        ==>     |         ==>         |
              .------ y         .------x'         .------x'      Y------.
             /                                   /                     /
          x /                                z''/                    Z/

         We thus see that the length along the second original direction
         (y) doesn't change stays in the second direction (Y). But the
         original first axis direction (x) is now represented by (Z). Note
         that we are talking about the semi-axis lengths (a scalar, not a
         vector), so +- is irrelevant. */
      else
        {
          naper[0] = oaper[2];
          naper[1] = oaper[1] / oaper[2];
          naper[2] = oaper[0] / oaper[2];
          naper[3] = naper[4] = naper[5] = 90;
        }

      break;

    case 6:
      if(oaper[1]>1 || oaper[2]>1)
        error(EXIT_FAILURE, 0, "atleast one of the second or third values "
              "to '--aperture' is larger than one. When size numbers are "
              "given to this option (in 3D), the second and third are the "
              "axis ratios (which must always be less than 1).");
      break;

    default:
      error(EXIT_FAILURE, 0, "%zu values given to '--aperture'. In 3D, "
            "this option can only take 1, 3, or 6 values. See the "
            "description of this option in the book for more with this "
            "command:\n\n"
            "    $ info %s\n", p->aperture->size, PROGRAM_EXEC);
    }

  /* If a new aperture was defined, then replace it with the exitsting
     one. */
  if(newaper)
    {
      gal_data_free(p->aperture);
      p->aperture=newaper;
    }
}





static size_t
ui_set_columns_sanity_check_read_aperture(struct matchparams *p)
{
  size_t ccol1n=0, ccol2n=0;

  /* Make sure the columns to match are given. */
  if(p->coord || p->kdtreemode==MATCH_KDTREE_BUILD)
    {
      if(p->ccol1==NULL)
        error(EXIT_FAILURE, 0, "no value given to '--ccol1' (necessary "
              "with '--coord')");
    }
  else
    {
      if(p->ccol1==NULL || p->ccol2==NULL)
        error(EXIT_FAILURE, 0, "both '--ccol1' and '--ccol2' must be "
              "given. They specify the columns containing the "
              "coordinates to match");
    }

  /* Make sure the same number of columns is given to both. Note that a
     second catalog is only necessary when we aren't building a k-d
     tree. */
  ccol1n = p->ccol1->size;
  if( p->kdtreemode!=MATCH_KDTREE_BUILD )
    {
      ccol2n = p->coord ? p->coord->size : p->ccol2->size;
      if(ccol1n!=ccol2n)
        error(EXIT_FAILURE, 0, "number of coordinates given to '--ccol1' "
              "(%zu) and '--%s' (%zu) must be equal.\n\n"
              "If you didn't call these options, run with '--checkconfig' "
              "to see which configuration file is responsible. You can "
              "always override the configuration file values by calling "
              "the option manually on the command-line", ccol1n,
              p->coord ? "coord" : "ccol2", ccol2n);
    }

  /* Read/check the aperture values. */
  if(p->aperture)
    switch(ccol1n)
      {
      case 1:
        if(p->aperture->size>1)
          error(EXIT_FAILURE, 0, "%zu values given to '--aperture'. In "
                "a 1D match, this option can only take one value",
                p->aperture->size);
        break;

      case 2: ui_read_columns_aperture_2d(p); break;
      case 3: ui_read_columns_aperture_3d(p); break;
      default:
        error(EXIT_FAILURE, 0, "%zu dimensional matches are not "
              "currently supported (maximum is 2 dimensions). The number "
              "of dimensions is deduced from the number of values given "
              "to '--ccol1' (or '--coord') and '--ccol2'", ccol1n);
      }
  else if ( p->kdtreemode != MATCH_KDTREE_BUILD )
    error(EXIT_FAILURE, 0, "no matching aperture specified. Please "
          "use the '--aperture' option to define the acceptable "
          "aperture for matching the coordinates (in the same units "
          "as each dimension). An aperture is not necessary when you "
          "are building a k-d tree or when you want an outer match; "
          "but this run is not any of those. Please run 'info %s' for "
          "more on the format of '--aperture'", PROGRAM_EXEC);

  /* Return the number of dimensions. */
  return ccol1n;
}





/* Save the manually given coordinates (with '--coord') in column format (a
   list of datasets). */
static gal_data_t *
ui_set_columns_from_coord(struct matchparams *p)
{
  size_t i, one=1;
  gal_data_t *out=NULL;
  double *coord=p->coord->array;

  /* Write each given value as a one-row table column (single element
     datasets that are linked) */
  for(i=0;i<p->coord->size;++i)
    {
      gal_list_data_add_alloc(&out, NULL, GAL_TYPE_FLOAT64, 1, &one, NULL,
                              0, -1, 1, NULL, NULL, NULL);
      *((double *)(out->array))=coord[i];
    }
  gal_list_data_reverse(&out);

  /* Return the list. */
  return out;
}





/* We want to keep the columns as double type. So what-ever their original
   type is, convert it. */
static gal_data_t *
ui_read_columns_to_double(struct matchparams *p, char *filename,
                          char *hdu, gal_list_str_t *cols,
                          size_t numcols, char *hdu_option_name)
{
  gal_data_t *tmp, *ttmp, *tout, *out=NULL;
  struct gal_options_common_params *cp=&p->cp;
  char *diff_cols_error="%s: the number of columns matched (%zu) "
    "differs from the number of usable calls to '--ccol1' (%zu). "
    "Please give more specific values to '--ccol1' (column "
    "numberes are the only identifiers guaranteed to be unique).";

  /* Read the columns. Note that the first input's name can be NULL (if the
     user intended to use the standrad input). Also note that this function
     is called more than one time, so if the Standard input is already read
     once, we don't want to write a blank list over it (the Standard input
     will be empty after being read). */
  if(p->stdinlines==NULL)
    p->stdinlines=gal_options_check_stdin(filename, p->cp.stdintimeout,
                                          "input");
  tout=gal_table_read(filename, hdu, filename ? NULL : p->stdinlines,
                      cols, cp->searchin, cp->ignorecase, cp->numthreads,
                      cp->minmapsize, p->cp.quietmmap, NULL,
                      hdu_option_name);

  /* A small sanity check. */
  if(gal_list_data_number(tout)!=numcols)
    error(EXIT_FAILURE, 0, diff_cols_error,
          gal_checkset_dataset_name(filename, hdu),
          gal_list_data_number(tout), numcols);

  /* Go over the columns and see if they are double or not. To keep things
     simple, we'll keep a new list even if all the types are float64.*/
  tmp=tout;
  while(tmp!=NULL)
    {
      /* If the coordinate column has more than one dimension (a vector
         column), the abort with an error: */
      if(tmp->ndim!=1)
        error(EXIT_FAILURE, 0, "only single-valued columns can be given "
              "to '--ccol1' or '--ccol2'");

      /* We need ot set the 'next' pointer  */
      ttmp=tmp->next;
      tmp->next=NULL;

      /* Correct the type if necessary. */
      if(tmp->type==GAL_TYPE_FLOAT64)
        gal_list_data_add(&out, tmp);
      else
        gal_list_data_add(&out,
                          gal_data_copy_to_new_type_free(tmp,
                                                         GAL_TYPE_FLOAT64) );

      /* Set 'tmp' to the initial 'next pointer. */
      tmp=ttmp;
    }

  /* The 'out' above is in reverse, so correct it and return */
  gal_list_data_reverse(&out);
  return out;
}




#define UI_READ_KDTREE_FIXED_MSG "You can construct a k-d tree "        \
  "for your first input table using the command below (just set your "  \
  "desired output name, and give that file to '--kdtree', in a later "  \
  "call with the same first input catalog):\n\n"                        \
  "    astmatch %s -h%s --ccol1=%s,%s --kdtree=build "                  \
  "--output=KDTREE.fits\n\n"                                            \
  "To learn more about the expected format of k-d trees in Gnuastro, "  \
  "please run this command: 'info gnuastro k-d'"

static void
ui_read_kdtree(struct matchparams *p)
{
  size_t *sizetarr;
  char **strarr1 = p->ccol1->array;
  gal_data_t *keysll=gal_data_array_calloc(1);

  /* Read the external k-d tree. */
  p->kdtreedata=gal_table_read(p->kdtree, p->kdtreehdu, NULL,
                               NULL, GAL_TABLE_SEARCH_NAME, 1,
                               p->cp.numthreads, p->cp.minmapsize,
                               p->cp.quietmmap, NULL, "--kdtreehdu");

  /* It has to have only two columns, with an unsigned 32-bit integer
     type in each. */
  if( gal_list_data_number(p->kdtreedata)!=2 )
    error(EXIT_FAILURE, 0, "%s (hdu: %s, that was given to "
          "'--kdtree') has %zu column(s), but should contain 2 "
          "columns. " UI_READ_KDTREE_FIXED_MSG,
          p->kdtree, p->kdtreehdu, gal_list_data_number(p->kdtreedata),
          p->input1name, p->cp.hdu, strarr1[0], strarr1[1]);
  if( p->kdtreedata->type!=GAL_TYPE_UINT32
      || p->kdtreedata->next->type!=GAL_TYPE_UINT32 )
    error(EXIT_FAILURE, 0, "%s (hdu: %s, that was given to "
          "'--kdtree') doesn't have unsigned 32-bit integer columns. "
          UI_READ_KDTREE_FIXED_MSG, p->kdtree, p->kdtreehdu,
          p->input1name, p->cp.hdu, strarr1[0], strarr1[1]);
  if( p->kdtreedata->size!=p->cols1->size )
    error(EXIT_FAILURE, 0, "%s (hdu: %s, that was given to "
          "'--kdtree') doesn't have the same number of rows as "
          "%s, which is the first given input. "
          UI_READ_KDTREE_FIXED_MSG, p->kdtree, p->kdtreehdu,
          gal_fits_name_save_as_string(p->input1name, p->cp.hdu),
          p->input1name, p->cp.hdu, strarr1[0], strarr1[1]);

  /* Read the k-d tree root. */
  keysll[0].type=GAL_TYPE_SIZE_T;
  keysll[0].name=MATCH_KDTREE_ROOT_KEY;
  gal_fits_key_read(p->kdtree, p->kdtreehdu, keysll, 0, 0, "--kdtreehdu");
  if(keysll[0].status)
    error(EXIT_FAILURE, 0, "%s (hdu: %s, that was given to "
          "'--kdtree') doesn't have the '%s' keyword, or it "
          "couldn't be read as a number. " UI_READ_KDTREE_FIXED_MSG,
          p->kdtree, p->kdtreehdu, MATCH_KDTREE_ROOT_KEY,
          p->input1name, p->cp.hdu, strarr1[0], strarr1[1]);
  sizetarr=keysll[0].array;
  p->kdtreeroot=sizetarr[0];

  /* Clean up: since the 'name' component wasn't allocated in this example,
     we should set it to NULL before calling 'gal_data_array_free'. */
  keysll[0].name=NULL;
  gal_data_array_free(keysll, 1, 1);
}





/* Read catalog columns */
static void
ui_read_columns(struct matchparams *p)
{
  size_t i, ndim;
  char **strarr1, **strarr2;
  gal_list_str_t *cols1=NULL, *cols2=NULL;

  /* Basic sanity checks and reading of aperture values. */
  ndim=ui_set_columns_sanity_check_read_aperture(p);

  /* Convert the array of strings to a list of strings for the column
     names. */
  strarr1 = p->ccol1->array;
  strarr2 = ( (p->coord || p->kdtreemode==MATCH_KDTREE_BUILD)
              ? NULL
              : p->ccol2->array );
  for(i=0;i<ndim;++i)
    {
      gal_list_str_add(&cols1, strarr1[i], 1);
      if(strarr2) gal_list_str_add(&cols2, strarr2[i], 1);
    }
  gal_list_str_reverse(&cols1);
  if(cols2) gal_list_str_reverse(&cols2);

  /* Read-in the columns. */
  p->cols1=ui_read_columns_to_double(p, p->input1name, p->cp.hdu,
                                     cols1, ndim, "--hdu");
  if( p->kdtreemode!=MATCH_KDTREE_BUILD )
    p->cols2=( p->coord
               ? ui_set_columns_from_coord(p)
               : ui_read_columns_to_double(p, p->input2name, p->hdu2,
                                           cols2, ndim, "--hdu2") );

  /* If an external k-d tree is given, read it and make sure it has the
     same number of rows as the first input and the proper datatype. */
  if( p->kdtreemode==MATCH_KDTREE_FILE )
    ui_read_kdtree(p);

  /* If we are in k-d tree based matching and the second dataset is smaller
     than the first, and we are not in an outer-match (where the order of
     the inputs matter) print a warning to let the user know that the speed
     can be greatly improved if they swap the two. */
  if( !p->cp.quiet
      && p->kdtreemode!=MATCH_KDTREE_BUILD
      && p->kdtreemode!=MATCH_KDTREE_DISABLE
      && p->cols1->size > (2*p->cols2->size)
      && p->type!=GAL_MATCH_ARRANGE_OUTER)
    error(EXIT_SUCCESS, 0, "TIP: the matching speed will improve "
          "if you swap the two inputs. Recall that for an inner or "
          "full match, the order of inputs does not matter. "
          "In the k-d tree based matching, a k-d tree will be "
          "constructed for the first input and parallelized matching "
          "will occur on the second input's rows. Fewer first-input "
          "rows therefore mean a smaller k-d tree (which is easier "
          "to build), and more second-input rows will be better "
          "distributed in your system's processing threads");

  /* Free the extra spaces. */
  gal_list_str_free(cols1, 1);
  gal_list_str_free(cols2, 1);
}





static void
ui_preparations_out_cols_read(gal_list_str_t **list, char *string)
{
  /* Make sure the string is not empty. */
  if(string[0]=='\0')
    error(EXIT_FAILURE, 0, "the value(s) given to '--outcols' should "
          "a column identifier (counter, name or '_all') after the 'a' "
          "or 'b'");

  /* Add the string to the list of columns for this column. */
  gal_list_str_add(list, string, 0);
}





static void
ui_preparations_out_cols(struct matchparams *p)
{
  void *rptr;
  int goodvalue;
  gal_data_t *read;
  uint8_t readtype;
  char *col, **strarr=p->outcols->array;
  size_t i, one=1, ndim=p->coord?p->coord->size:0;

  /* Go over all the values and put the respective column identifier in the
     proper list. */
  for(i=0;i<p->outcols->size;++i)
    {
      /* For easy reading. */
      col=strarr[i];

      /* In no-match mode, the same column will be used from both catalogs
         so things are easier. */
      if(p->notmatched)
        {
          gal_list_str_add(&p->acols, col, 0);
          gal_list_str_add(&p->bcols, col, 0);
        }

      /* In match mode, we need to know which column should come from which
         catalog. */
      else
        switch(col[0])
          {
          case 'a':
            ui_preparations_out_cols_read(&p->acols, col+1);
            break;
          case 'b':
            /* With '--coord', only numbers that are smaller than the
               number of the dimensions are acceptable. */
            if(p->coord)
              {
                goodvalue=0;
                rptr=gal_type_string_to_number(col+1, &readtype);
                if(rptr)
                  {
                    read=gal_data_alloc(rptr, readtype, 1, &one, NULL, 0,
                                        -1, 1, NULL, NULL, NULL);
                    if(gal_type_is_int(readtype))
                      {
                        read=gal_data_copy_to_new_type_free(read,
                                                            GAL_TYPE_LONG);
                        if( *((long *)(read->array)) <= ndim )
                          goodvalue=1;
                      }
                    gal_data_free(read);
                  }
                if(goodvalue==0)
                  error(EXIT_FAILURE, 0, "bad value to second catalog "
                        "column (%s) of '--outcols'.\n\n"
                        "With the '--coord' option, the second catalog "
                        "is assumed to have a single row and the given "
                        "number of columns. Therefore when using "
                        "'--outcols', only integers that are less than "
                        "the number of dimensions (%zu in this case) "
                        "are acceptable", col+1, ndim);
              }
            ui_preparations_out_cols_read(&p->bcols, col+1);
            break;
          default:
            error(EXIT_FAILURE, 0, "'%s' is not a valid value for "
                  "'--outcols'.\n\n"
                  "The first character of each value to this option "
                  "must be either 'a' or 'b'. The former specifies a "
                  "column from the first input and the latter a "
                  "column from the second. The characters after them "
                  "can be any column identifier (number, name, or "
                  "regular expression). For more on column selection, "
                  "please run this command:\n\n"
                  "    $ info gnuastro \"Selecting table columns\"\n",
                  col);
          }
    }

  /* Reverse the lists so they correspond to the input order. */
  gal_list_str_reverse(&p->acols);
  gal_list_str_reverse(&p->bcols);
}





#define UI_OUT_SUFFNAME "-matched"
static void
ui_preparations_out_name(struct matchparams *p)
{
  /* To temporarily keep the original value. */
  char *suffix;
  uint8_t keepinputdir_orig;
  char *refname = p->input1name ? p->input1name : p->input2name;

  /* Set the output file(s) name(s). */
  if(p->logasoutput)
    {
      if(p->cp.output)
        gal_checkset_allocate_copy(p->cp.output, &p->logname);
      else
        {
          if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
            p->logname=gal_checkset_automatic_output(&p->cp, refname,
                                                     UI_OUT_SUFFNAME".txt");
          else
            p->logname=gal_checkset_automatic_output(&p->cp, refname,
                                                     UI_OUT_SUFFNAME".fits");
        }

      /* Make sure a file with this name doesn't exist. */
      gal_checkset_writable_remove(p->logname, NULL, 0, p->cp.dontdelete);

      /* The main output name needs to be available in p->out1name (for the
         final step when we want to write the input configurations as FITS
         keywords). */
      gal_checkset_allocate_copy(p->logname, &p->out1name);
    }
  else
    {
      if(p->outcols || p->coord || p->kdtreemode==MATCH_KDTREE_BUILD)
        {
          if(p->cp.output)
            gal_checkset_allocate_copy(p->cp.output, &p->out1name);
          else
            {
              suffix = ( p->kdtreemode==MATCH_KDTREE_BUILD
                         ? "-kdtree.fits"
                         : ( p->cp.tableformat==GAL_TABLE_FORMAT_TXT
                             ? UI_OUT_SUFFNAME".txt"
                             : UI_OUT_SUFFNAME".fits") );
              p->out1name = gal_checkset_automatic_output(&p->cp,
                                                          refname, suffix);
            }
          gal_checkset_writable_remove(p->out1name, refname, 0,
                                       p->cp.dontdelete);
        }
      else
        {
          /* Set 'p->out1name' and 'p->out2name'. */
          if(p->cp.output)
            {
              if( gal_fits_name_is_fits(p->cp.output) )
                {
                  gal_checkset_allocate_copy(p->cp.output, &p->out1name);
                  gal_checkset_allocate_copy(p->cp.output, &p->out2name);
                }
              else
                {
                  /* Here, we are be using the output name as input to the
                     automatic output generating function (usually it is
                     the input name, not the output name). Therefore, the
                     'keepinputdir' variable should be 1. So we will
                     temporarily change it here, then set it back to what
                     it was. */
                  keepinputdir_orig=p->cp.keepinputdir;
                  p->cp.keepinputdir=1;
                  p->out1name=gal_checkset_automatic_output(&p->cp,
                                                            p->cp.output,
                                                            UI_OUT_SUFFNAME
                                                            "-1.txt");
                  p->out2name=gal_checkset_automatic_output(&p->cp,
                                                            p->cp.output,
                                                            UI_OUT_SUFFNAME
                                                            "-2.txt");
                  p->cp.keepinputdir=keepinputdir_orig;
                }
            }
          else
            {
              if(p->cp.tableformat==GAL_TABLE_FORMAT_TXT)
                {
                  p->out1name=gal_checkset_automatic_output(&p->cp, refname,
                                                            UI_OUT_SUFFNAME
                                                            "-1.txt");
                  p->out2name=gal_checkset_automatic_output(&p->cp,
                                                            p->input2name,
                                                            UI_OUT_SUFFNAME
                                                            "-2.txt");
                }
              else
                {
                  p->out1name=gal_checkset_automatic_output(&p->cp, refname,
                                                            UI_OUT_SUFFNAME
                                                            ".fits");
                  gal_checkset_allocate_copy(p->out1name, &p->out2name);
                }
            }

          /* Make sure no file with these names exists. */
          gal_checkset_writable_remove(p->out1name, refname, 0,
                                       p->cp.dontdelete);
          gal_checkset_writable_remove(p->out2name, refname, 0,
                                       p->cp.dontdelete);
        }

      /* If a log file is necessary, set its name here. */
      if(p->cp.log)
        p->logname=gal_checkset_automatic_output(&p->cp, p->out1name,
                                                 UI_OUT_SUFFNAME
                                                 "-log.fits");
    }
}





static void
ui_preparations(struct matchparams *p)
{
  int intype;

  /* Since Match only works on tables, 'searchin' is necessary. */
  if(p->cp.searchin==0)
    error(EXIT_FAILURE, 0, "no '--searchin' option specified. Please "
          "run the following command for more information:\n\n"
          "    $ info gnuastro \"selecting table columns\"\n");

  /* We will base the mode on the first input, then check with the
     second. Note that when the first is from standard input (it is
     'NULL'), then we go into catalog mode because currently we assume
     standard input is only for plain text and WCS matching is not defined
     on plain text. */
  if( p->input1name && gal_fits_file_recognized(p->input1name) )
    {
      intype=gal_fits_hdu_format(p->input1name, p->cp.hdu, "--hdu");
      if(intype==IMAGE_HDU)
        error(EXIT_FAILURE, 0, "%s: not a catalog",
              gal_checkset_dataset_name(p->input1name, p->cp.hdu));
    }

  /* Same check for second input. */
  if(p->input2name && gal_fits_file_recognized(p->input2name))
    {
      intype=gal_fits_hdu_format(p->input2name, p->hdu2, "--hdu2");
      if(intype==IMAGE_HDU)
        error(EXIT_FAILURE, 0, "%s: not a catalog",
              gal_checkset_dataset_name(p->input2name, p->hdu2));
    }

  /* Read the input columns and prepare output columns. */
  ui_read_columns(p);
  if(p->outcols)
    ui_preparations_out_cols(p);

  /* Set the output filename. */
  ui_preparations_out_name(p);
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
ui_read_check_inputs_setup(int argc, char *argv[], struct matchparams *p)
{
  size_t nthreads;
  struct gal_options_common_params *cp=&p->cp;


  /* Include the parameters necessary for argp from this program ('args.h')
     and for the common options to all Gnuastro ('commonopts.h'). We want
     to directly put the pointers to the fields in 'p' and 'cp', so we are
     simply including the header here to not have to use long macros in
     those headers which make them hard to read and modify. This also helps
     in having a clean environment: everything in those headers is only
     available within the scope of this function. */
#include <gnuastro-internal/commonopts.h>
#include "args.h"


  /* Initialize the options and necessary information.  */
  ui_initialize_options(p, program_options, gal_commonopts_options);


  /* Read the command-line options and arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");


  /* Read the configuration files and set the common values. */
  gal_options_read_config_set(&p->cp);


  /* Sanity check only on options. */
  ui_check_only_options(p);


  /* Print the option values if asked. Note that this needs to be done
     after the option checks so un-sane values are not printed in the
     output state. */
  gal_options_print_state(&p->cp);


  /* Check that the options and arguments fit well with each other. Note
     that arguments don't go in a configuration file. So this test should
     be done after (possibly) printing the option values. */
  ui_check_options_and_arguments(p);


  /* Read/allocate all the necessary starting arrays. */
  ui_preparations(p);


  /* If the output is a FITS table, prepare all the options as FITS
     keywords to write in output later. */
  if(gal_fits_name_is_fits(p->out1name))
    gal_options_as_fits_keywords(&p->cp);


  /* Let the user know that processing has started. */
  if(!p->cp.quiet)
    {
      printf(PROGRAM_NAME" "PACKAGE_VERSION" started on %s",
             ctime(&p->rawtime));
      nthreads=p->kdtreemode==MATCH_KDTREE_DISABLE ? 1 : p->cp.numthreads;
      printf("  - Using %zu CPU thread%s%s\n", nthreads,
             nthreads==1 ? "." : "s.",
             ( p->kdtreemode==MATCH_KDTREE_DISABLE
               ? " (sort-based match only uses a single thread)" : ""));
      printf("  - Match algorithm: %s\n",
             p->kdtree ? "k-d tree" : "sort-based");
      printf("  - Input-1: %s; %zu rows\n",
             gal_fits_name_save_as_string(p->input1name, p->cp.hdu),
             p->cols1->size);
      if(p->kdtreemode==MATCH_KDTREE_FILE)
        printf("  - Input-1 k-d tree: %s\n",
               gal_fits_name_save_as_string(p->kdtree, p->kdtreehdu));
      if(p->kdtreemode!=MATCH_KDTREE_BUILD)
        printf("  - Input-2: %s; %zu rows\n",
               p->coord ? "from --coord"
               : gal_fits_name_save_as_string(p->input2name, p->hdu2),
               p->cols2->size);
    }
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
ui_free_report(struct matchparams *p, struct timeval *t1)
{
  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(p->aperture);
  free(p->out1name);
  free(p->out2name);
  free(p->cp.output);
  gal_data_free(p->ccol1);
  gal_data_free(p->ccol2);
  gal_data_free(p->outcols);
  gal_list_data_free(p->cols1);
  gal_list_data_free(p->cols2);
  if(p->kdtree) free(p->kdtree);
  gal_list_str_free(p->acols, 0);
  gal_list_str_free(p->bcols, 0);
  if(p->kdtreehdu) free(p->kdtreehdu);
  gal_list_str_free(p->stdinlines, 1);

  /* Print the final message. */
  if(!p->cp.quiet)
    gal_timing_report(t1, PROGRAM_NAME" finished in: ", 0);
}
