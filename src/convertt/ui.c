/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (gnuastro) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <math.h>
#include <argp.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>

#include "timing.h"	/* Includes time.h and sys/time.h   */
#include "txtarrayvv.h"
#include "linkedlist.h"
#include "configfiles.h"
#include "fitsarrayvv.h"

#include "main.h"

#include "eps.h"
#include "jpeg.h"
#include "args.h"



/* Set the file names of the places where the default parameters are
   put. */
#define CONFIG_FILE SPACK CONF_POSTFIX
#define SYSCONFIG_FILE SYSCONFIG_DIR CONFIG_FILE
#define USERCONFIG_FILEEND USERCONFIG_DIR CONFIG_FILE
#define CURDIRCONFIG_FILE CURDIRCONFIG_DIR CONFIG_FILE





/**************************************************************/
/**************       Options and parameters    ***************/
/**************************************************************/
void
readconfig(char *filename, struct converttparams *p)
{
  int tmp;
  FILE *fp;
  size_t lineno=0, len=200;
  char *line, *name, *value;
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;
  char key='a';	/* Not used, just a place holder. */

  /* When the file doesn't exist or can't be opened, it is ignored. It
     might be intentional, so there is no error. If a parameter is
     missing, it will be reported after all defaults are read. */
  fp=fopen(filename, "r");
  if (fp==NULL) return;


  /* Allocate some space for `line` with `len` elements so it can
     easily be freed later on. The value of `len` is arbitarary at
     this point, during the run, getline will change it along with the
     pointer to line. */
  errno=0;
  line=malloc(len*sizeof *line);
  if(line==NULL)
    error(EXIT_FAILURE, errno, "ui.c: %lu bytes in readdefaults",
	  len * sizeof *line);

  /* Read the tokens in the file:  */
  while(getline(&line, &len, fp) != -1)
    {
      /* Prepare the "name" and "value" strings, also set lineno. */
      STARTREADINGLINE;

      /* Inputs: */
      if(strcmp(name, "hdu")==0)
	{
	  if(cp->hduset) continue;
	  errno=0;
	  cp->hdu=malloc(strlen(value)+1);
	  if(cp->hdu==NULL)
	    error(EXIT_FAILURE, errno, "Space for HDU.");
	  strcpy(cp->hdu, value);
	  cp->hduset=1;
	}
      else if(strcmp(name, "hdu2")==0)
	{
	  if(up->hdu2set) continue;
	  errno=0;
	  up->hdu2=malloc(strlen(value)+1);
	  if(up->hdu2==NULL)
	    error(EXIT_FAILURE, errno, "Space for hdu2.");
	  strcpy(up->hdu2, value);
	  up->hdu2set=1;
	}
      else if(strcmp(name, "hdu3")==0)
	{
	  if(up->hdu3set) continue;
	  errno=0;
	  up->hdu3=malloc(strlen(value)+1);
	  if(up->hdu3==NULL)
	    error(EXIT_FAILURE, errno, "Space for hdu3.");
	  strcpy(up->hdu3, value);
	  up->hdu3set=1;
	}
      else if(strcmp(name, "hdu4")==0)
	{
	  if(up->hdu4set) continue;
	  errno=0;
	  up->hdu4=malloc(strlen(value)+1);
	  if(up->hdu4==NULL)
	    error(EXIT_FAILURE, errno, "Space for hdu4.");
	  strcpy(up->hdu4, value);
	  up->hdu4set=1;
	}





      /* Outputs: */
      else if(strcmp(name, "output")==0)
	{
	  if(cp->outputset) continue;
	  errno=0;
	  cp->output=malloc(strlen(value)+1);
	  if(cp->output==NULL)
	    error(EXIT_FAILURE, errno, "Space for output");
	  strcpy(cp->output, value);
	  cp->outputset=1;
	}
      else if(strcmp(name, "quality")==0)
	{
	  if(up->qualityset) continue;
          intsmallerequalto(value, &p->quality, name, key,
                            p->cp.spack, filename, lineno, 100);
          if(p->quality<0)
            error(EXIT_FAILURE, 0, "The quality option should be positive.");
	  up->qualityset=1;
	}
      else if(strcmp(name, "widthincm")==0)
	{
	  if(up->widthincmset) continue;
          floatl0(value, &p->widthincm, name, key, SPACK, filename, lineno);
	  up->widthincmset=1;
	}
      else if(strcmp(name, "borderwidth")==0)
	{
	  if(up->borderwidthset) continue;
          intelzero(value, &p->borderwidth, name, key, SPACK,
                    filename, lineno);
	  up->borderwidthset=1;
	}





      /* Flux: */
      else if(strcmp(name, "fluxlow")==0)
	{
	  if(up->fluxlowset) continue;
          anydouble(value, &p->fluxlow, name, key, p->cp.spack,
                   filename, lineno);
          up->fluxlowset=1;
	}
      else if(strcmp(name, "fluxhigh")==0)
	{
	  if(up->fluxhighset) continue;
          anydouble(value, &p->fluxhigh, name, key, p->cp.spack,
                   filename, lineno);
          up->fluxhighset=1;
	}
      else if(strcmp(name, "maxbyte")==0)
	{
	  if(up->maxbyteset) continue;
          intsmallerequalto(value, &tmp, "maxbyte", key,
                            p->cp.spack, NULL, 0, UINT8_MAX);
          if(tmp<0)
            error(EXIT_FAILURE, 0, "--maxbyte (-m) should be positive.");
          p->maxbyte=tmp;
          p->up.maxbyteset=1;
	}



      /* Operating modes: */
      else if(strcmp(name, "numthreads")==0)
	{
	  if(cp->numthreadsset) continue;
	  sizetlzero(value, &cp->numthreads, name, key, SPACK,
		     filename, lineno);
	  cp->numthreadsset=1;
	}


      else
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      "`%s` not recognized.\n", name);
    }

  free(line);
  fclose(fp);
}





void
printvalues(FILE *fp, struct converttparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  /* Print all the options that are set. Separate each group with a
     commented line explaining the options in that group. */
  fprintf(fp, "\n# Input image:\n");
  if(cp->hduset)
    {
      if(stringhasspace(cp->hdu))
	fprintf(fp, CONF_SHOWFMT"\"%s\"\n", "hdu", cp->hdu);
      else
	fprintf(fp, CONF_SHOWFMT"%s\n", "hdu", cp->hdu);
    }
  if(up->hdu2set)
    {
      if(stringhasspace(up->hdu2))
	fprintf(fp, CONF_SHOWFMT"\"%s\"\n", "hdu2", up->hdu2);
      else
	fprintf(fp, CONF_SHOWFMT"%s\n", "hdu2", up->hdu2);
    }
  if(up->hdu3set)
    {
      if(stringhasspace(up->hdu3))
	fprintf(fp, CONF_SHOWFMT"\"%s\"\n", "hdu3", up->hdu3);
      else
	fprintf(fp, CONF_SHOWFMT"%s\n", "hdu3", up->hdu3);
    }
  if(up->hdu4set)
    {
      if(stringhasspace(up->hdu4))
	fprintf(fp, CONF_SHOWFMT"\"%s\"\n", "hdu4", up->hdu4);
      else
	fprintf(fp, CONF_SHOWFMT"%s\n", "hdu4", up->hdu4);
    }


  fprintf(fp, "\n# Output parameters:\n");
  if(cp->outputset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "output", cp->output);
  if(up->qualityset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "quality", p->quality);
  if(up->widthincmset)
    fprintf(fp, CONF_SHOWFMT"%.2f\n", "widthincm", p->widthincm);
  if(up->borderwidthset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "borderwidth", p->borderwidth);


  fprintf(fp, "\n# Output flux display:\n");
  if(up->fluxlowset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "fluxlow", p->fluxlow);
  if(up->fluxhighset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "fluxhigh", p->fluxhigh);
  if(up->maxbyteset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "maxbyte", (int)(p->maxbyte));
}





/* Note that numthreads will be used automatically based on the
   configure time. */
void
checkifset(struct converttparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  int intro=0;
  if(cp->hduset==0)
    REPORT_NOTSET("hdu");
  if(up->hdu2set==0)
    REPORT_NOTSET("hdu2");
  if(up->hdu3set==0)
    REPORT_NOTSET("hdu3");
  if(up->hdu4set==0)
    REPORT_NOTSET("hdu4");
  if(cp->outputset==0)
    REPORT_NOTSET("output");
  if(up->qualityset==0)
    REPORT_NOTSET("quality");
  if(up->widthincmset==0)
    REPORT_NOTSET("widthincm");
  if(up->borderwidthset==0)
    REPORT_NOTSET("borderwidth");
  if(up->fluxlowset==0)
    REPORT_NOTSET("fluxlow");
  if(up->fluxhighset==0)
    REPORT_NOTSET("fluxhigh");
  if(up->maxbyteset==0)
    REPORT_NOTSET("maxbyte");

  END_OF_NOTSET_REPORT;
}



















/****************************************************************
 *****************   Read convert values:    ********************
 ****************************************************************/
struct change *
makechangestruct(char *arg)
{
  char *p=arg;
  struct change *out=NULL, *c;

  while(*p!='\0')
    {
      while(*p==' ') {++p; continue;} /* Skip all space characters. */
      errno=0;
      c=malloc(sizeof *c);
      if(c==NULL) error(EXIT_FAILURE, 0, "%lu bytes for struct change",
                        sizeof *c);
      c->from=strtof(p, &p);
      while(*p==' ') {++p; continue;}
      if(*p==':') ++p;
      else
	{
	  fprintf(stderr, PACKAGE": In the conversion option, [from] "
		  "and [to] values should be separated by a ':'. You "
		  "have given a '%c': %s\n", *p, arg);
	  exit(EXIT_FAILURE);
	}
      c->to=strtof(p, &p);
      while(*p==' ') {++p; continue;}
      if(*p==',') p++;
      else if(*p!='\0')
	{
	  fprintf(stderr, PACKAGE": In the conversion option, [from] "
		  "and [to] pairs should be separated by a ','. You have "
		  "provided a '%c': %s\n", *p, arg);
	  exit(EXIT_FAILURE);
	}
      c->next=out;
      out=c;
    }
  /*
  {
    struct change *tmp;
    for(tmp=out;tmp!=NULL;tmp=tmp->next)
      printf("%f --> %f\n", tmp->from, tmp->to);
  }
  */
  return out;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* We know cp->output is a known suffix, we just don't know if it has
   a `.` before it or not. If it doesn't one will be added to it an
   the output name will be set. */
void
adddotautomaticoutput(struct converttparams *p)
{
  size_t i;
  struct commonparams *cp=&p->cp;
  char *tmp, *basename="output.txt";

  /* Find the first file name in the input(s). */
  for(i=0;i<p->numinputs;++i)
    if(strcmp(p->names[i], BLANKCHANNELNAME))
      {
        basename=p->names[i];
        break;
      }

  /* If the suffix does not start with a `.`, put one there. */
  if(cp->output[0]!='.')
    {
      errno=0;
      tmp=malloc(strlen(cp->output)+10*sizeof *tmp);
      if(tmp==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for suffix name.",
              strlen(cp->output)+10*sizeof *tmp);
      sprintf(tmp, ".%s", cp->output);
      free(cp->output);
      cp->output=tmp;
    }

  /* Set the automatic output and make sure we have write access. */
  automaticoutput(basename, cp->output, cp->removedirinfo,
                  cp->dontdelete, &cp->output);
  if( dir0file1(cp->output, cp->dontdelete)==0 )
    error(EXIT_FAILURE, 0, "%s is a directory.", cp->output);
}





void
sanitycheck(struct converttparams *p)
{
  size_t i, j;
  struct commonparams *cp=&p->cp;

  /* The flux range: */
  if(p->fluxlow>p->fluxhigh)
    error(EXIT_FAILURE, 0, "The value of `--fluxlow` (`-L`, %.3f) is "
          "larger than `--fluxhigh` (`-H`, %.3f).", p->fluxlow, p->fluxhigh);

  /* Make sure there are 1 (for grayscale), 3 (for RGB) or 4 (for
     CMYK) color channels. */
  if(p->numch!=1 && p->numch!=3 && p->numch!=4)
    error(EXIT_FAILURE, 0, "The number of input color channels has to be "
          "1 (for non image data, grayscale or only K channel in CMYK), "
          "3 (for RGB) and 4 (for CMYK). You have given %lu color channels. "
          "Note that some file formats (for example JPEG) can contain more "
          "than one color channel.", p->numch);

  /* Make sure that there is atleast one input file (not only blank)
     and set the sizes of the blank channels to the first non-blank
     value. */
  for(i=0;i<p->numch;++i)
    if(p->isblank[i]==0) break;
  if(i==p->numch)
    error(EXIT_FAILURE, 0, "All the input(s) are of type blank!");
  for(j=0;j<p->numch;++j)
    if(p->isblank[j])
      {
        p->s0[j]=p->s0[i];
        p->s1[j]=p->s1[i];
      }

  /* Check if all the input sources have the same size: */
  if(p->numch>1)
    {
      for(i=1;i<p->numch;++i)
        if(p->s0[i]!=p->s0[0] || p->s1[i]!=p->s1[0])
          break;
      if(i!=p->numch)
        {
          for(i=0;i<p->numch;++i)
            fprintf(stderr, "Channel %lu is %lu x %lu pixels.\n", i,
                    p->s1[i], p->s0[i]);
          error(EXIT_FAILURE, 0, "The input color channels have different "
                "sizes.");
        }
    }

  /* Allocate space for the blank channels and set them to zero: */
  for(i=0;i<p->numch;++i)
    if(p->isblank[i])
      {
        errno=0;
        p->ch[i]=calloc(p->s0[0]*p->s1[0], sizeof *p->ch[i]);
        if(p->ch[i]==NULL)
          error(EXIT_FAILURE, errno, "Allocating %lu bytes for the blank "
                "channel %lu", p->s0[0]*p->s1[0]*sizeof *p->ch[i], i);
      }

  /* The output file name. First find the first non-blank file name: */
  if(nameisfits(cp->output))
    {
      p->outputtype=FITSFORMAT;
      if( nameisfitssuffix(cp->output) )
        adddotautomaticoutput(p);
    }
  else if(nameisjpeg(cp->output))
    {
      p->outputtype=JPEGFORMAT;
      if( nameisjpegsuffix(cp->output) )
        adddotautomaticoutput(p);
    }
  else if(nameiseps(cp->output))
    {
      p->outputtype=EPSFORMAT;
      if( nameisepssuffix(cp->output) )
        adddotautomaticoutput(p);
    }
  else if(nameispdf(cp->output))
    {
      p->outputtype=PDFFORMAT;
      if( nameispdfsuffix(cp->output) )
        adddotautomaticoutput(p);
    }
  else
    {
      /* If the length of the name is shorter than 4 characters, it is
         most probably a mis-spelled extension, warn the user. */
      if(strlen(cp->output)<=5)
        fprintf(stderr, SPACK": (Warning) Your output file name is `%s`, "
                "based on its length, it might be a mis-spelled extension. "
                "Your input is converted to a plain text format file with "
                "That name.",
                cp->output);

      p->outputtype=TXTFORMAT;

      /* If output type is not an image, there should only be one color
         channel: */
      if(p->numch>1)
        error(EXIT_FAILURE, 0, "Text output (`--output=%s`) can only be "
              "completed with one input color channel. You have given %lu. "
              "Note that some formats (for example JPEG) can have more than "
              "one color channel in each file. You can first convert the "
              "file to FITS, then convert the desired channel to text by "
              "specifying the HDU.",
              cp->output, p->numch);
    }

}




















/**************************************************************/
/***************        Preparations        *******************/
/**************************************************************/
void
preparearrays(struct converttparams *p)
{
  size_t i;
  void *array;
  double *d, *df;
  struct stll *tmp;
  char *hdu, **names=p->names;

  /* Put the names in the correct order. */
  i=p->numinputs-1;
  for(tmp=p->inputnames; tmp!=NULL; tmp=tmp->next)
    names[i--]=tmp->v;

  p->numch=0;
  for(i=0;i<p->numinputs;++i)
    {
      /* Check if p->numch has not exceeded 4. */
      if(p->numch>=4)
        error(EXIT_FAILURE, 0, "The number of input color channels (not "
              "files) has exceeded 4! Note that one file can contain more "
              "than one color channel.");

      /* FITS: */
      if( nameisfits(names[i]) )
        {
          switch(p->numch) /* Get the HDU value for this channel. */
            {
            case 0: hdu=p->cp.hdu; break;   case 1: hdu=p->up.hdu2; break;
            case 2: hdu=p->up.hdu3; break;  case 3: hdu=p->up.hdu4; break;
            default: error(EXIT_FAILURE, 0, "A bug! In parsing the input "
                           "FITS files, it has gone beyond four! Please "
                           "contact us so we can see what caused this "
                           "problem and fix it.");
            }
          p->numnul[p->numch]=fitsimgtoarray(names[i], hdu,
                                             &p->bitpixs[p->numch], &array,
                                             &p->s0[p->numch],
                                             &p->s1[p->numch]);
          changetype(array, p->bitpixs[p->numch],
                     p->s0[p->numch]*p->s1[p->numch], p->numnul[p->numch],
                     (void **)(&p->ch[p->numch]), DOUBLE_IMG);
          free(array);
          ++p->numch;
        }



      /* JPEG: */
      else if ( nameisjpeg(names[i]) )
        preparejpeg(p, names[i]);



      /* Blank: */
      else if(strcmp(names[i], BLANKCHANNELNAME)==0)
        {
          p->isblank[p->numch]=1;
          p->bitpixs[p->numch]=BYTE_IMG;
          ++p->numch;
        }



      /* EPS:  */
      else if ( nameiseps(names[i]) )
        error(EXIT_FAILURE, 0, "EPS files cannot be used as input. Since "
              "EPS files are not raster graphics, they are only used as "
              "output.");



      /* PDF:  */
      else if ( nameispdf(names[i]) )
        error(EXIT_FAILURE, 0, "PDF files cannot be used as input. Since "
              "PDF files are not raster graphics, they are only used as "
              "output.");


      /* Text: */
      else
        {
          txttoarray(names[i], &p->ch[p->numch],
                     &p->s0[p->numch], &p->s1[p->numch]);
          df = (d=p->ch[p->numch]) + p->s0[p->numch]*p->s1[p->numch];
          do if(isnan(*d++)) break; while(d<df);
          if(d==df)
            checkremovefile(TXTARRAYVVLOG, 0);
          else
            error(EXIT_FAILURE, 0, "%s contains non-numeric data, see %s.",
                  names[i], TXTARRAYVVLOG);
          p->bitpixs[p->numch]=DOUBLE_IMG;
          ++p->numch;
        }
    }
}


















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct converttparams *p)
{
  struct commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = DP_NUMTHREADS;
  cp->removedirinfo = 1;
  p->invert         = 1;

  /* Read the arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "Parsing arguments");

  /* Add the user default values and save them if asked. */
  CHECKSETCONFIG;

  /* Check if all the required parameters are set. */
  checkifset(p);

  /* Print the values for each parameter. */
  if(cp->printparams)
    REPORT_PARAMETERS_SET;

  /* Prepare the arrays: */
  preparearrays(p);

  /* Do a sanity check. */
  sanitycheck(p);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct converttparams *p)
{
  size_t i;
  struct stll *tmp, *ttmp;

  free(p->cp.hdu);
  free(p->up.hdu2);
  free(p->up.hdu3);
  free(p->up.hdu4);
  free(p->cp.output);
  for(i=0;i<4;++i) free(p->ch[i]);

  /* Free the input file names: */
  tmp=p->inputnames;
  while(tmp!=NULL)
    {
      ttmp=tmp->next;
      free(tmp);
      tmp=ttmp;
    }
}
