/*********************************************************************
NoiseChisel - Detect and segment signal in noise.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2015, Free Software Foundation, Inc.

Gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <stdlib.h>

#include "forqsort.h"
#include "neighbors.h"
#include "fitsarrayvv.h"

#include "main.h"

#include "label.h"
#include "clumps.h"
#include "segmentation.h"



/******************************************************************/
/*****************     Growing preparations     *******************/
/******************************************************************/
/* Set the growth threshold. We need one standard deviation for the
   full detection so we take the light weighted center of the full
   detected region and use that to find a standard deviation. Note
   that we are using the sky subtracted input image so gradients don't
   harm our results and that we only use pixels above the sky
   value. */
void
prepfirstgrowth(struct clumpsthreadparams *ctp)
{
  struct meshparams *smp=&ctp->p->smp;

  size_t *ind, *indf, is1=smp->s1;
  double x=0.0f, y=0.0f, fluxsum=0.0f;
  float growlimit, *imgss=ctp->p->imgss;
  long *olab=ctp->p->olab, *clab=ctp->p->clab;

  /* Try to find the flux weighted center (only on the pixels with
     positive flux) */
  indf=(ind=ctp->inds)+ctp->area;
  do
    if(imgss[*ind]>0.0f)
      {
        fluxsum+=imgss[*ind];
        x+=(*ind/is1)*imgss[*ind];
        y+=(*ind%is1)*imgss[*ind];
      }
  while(++ind<indf);

  /* Calculate the center, if no pixels were positive, use the
     geometric center (irrespective of flux). */
  if(fluxsum==0.0f)
    {
      ind=ctp->inds;
      do { x+=*ind/is1; y+=*ind%is1; } while(++ind<indf);
      x/=ctp->area;
      y/=ctp->area;
    }
  else
    {
      x/=fluxsum;
      y/=fluxsum;
    }

  /* First find the standard deviation on this detection, then use
     it to calculate the growth threshold. */
  ctp->std=smp->garray2[imgxytomeshid(smp, x, y)];
  growlimit=ctp->p->gthresh * ctp->std;

  /* Allocate and fill the array for the blank pixels and reset the
     olab pixels while you are at it. Note that after removing false
     clumps, usually most of the objects become blank, so we can just
     allocate an area the size of the detection: */
  errno=0;
  ctp->blankinds=malloc(ctp->area*sizeof *ctp->blankinds);
  if(ctp->blankinds==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for ctp->blankpixels in "
          "growclumps (clumps.c)", ctp->area*sizeof *ctp->blankinds);
  ctp->numblanks=0;
  indf=(ind=ctp->inds)+ctp->area;
  do
    {
      olab[*ind]=clab[*ind];
      if(olab[*ind]==SEGMENTINIT && imgss[*ind]>growlimit)
        ctp->blankinds[ctp->numblanks++]=*ind;
    }
  while(++ind<indf);
}




















/******************************************************************/
/*****************            Objects           *******************/
/******************************************************************/
void
thisdetectionisoneobject(struct clumpsthreadparams *ctp)
{
  size_t *ind, *indf;
  long *olab=ctp->p->olab;

  indf=(ind=ctp->inds)+ctp->area;
  do olab[*ind]=1; while(++ind<indf);
}





/* Find the adjacency matrixs (number, sum and signal to noise) for
   the rivers between potentially separate objects in a detection
   region. They have to be allocated prior to entering this funciton.

   The way to find connected objects is through an adjacency
   matrix. It is a square array with a side equal to numobjs. So to
   see if regions `a` and `b` are connected. All I have to do is to
   look at element a*numobjs+b or b*numobjs+a and get my answer. Since
   the number of objects in a given region will not be too high, this
   is efficient. */
void
adjacencymatrixs(struct clumpsthreadparams *ctp,
		 double *sns, double *sums, int *nums)
{
  int rpnum;
  float *imgss=ctp->p->imgss;
  size_t numclumps=ctp->numclumps;
  long wngb[WNGBSIZE], *olab=ctp->p->olab;
  size_t i, j, ii, *n, *nf, ngb[8], numngb;
  double ave, rpave, c=sqrt(1/ctp->p->cpscorr);
  size_t *ind, *indf, is0=ctp->p->lmp.s0, is1=ctp->p->lmp.s1;
  double err=ctp->std*ctp->std*(ctp->p->skysubtracted?1.0f:2.0f);

  /* Go over all the still-unlabeled pixels and see which labels they
     touch. */
  indf = (ind=ctp->inds) + ctp->area;
  do
    if( olab[*ind]==SEGMENTINIT && !isnan(imgss[*ind]) )
      {
	/* Initialize the values to be used: */
	ii=0;
	rpnum=1;
	rpave=imgss[*ind];
	memset(wngb, 0, sizeof(wngb));

	/* Find which grown clumps this river pixel touches.      */
	FILL_NGB_8_ALLIMG;
	nf=(n=ngb)+numngb;
	do
	  if( olab[*n]>0 )
	    {
              if(!isnan(imgss[*n])) /* Add the flux of this neighbor */
                {                   /* pixel for finding the average */
                  ++rpnum;          /* of this river pixel later.    */
                  rpave+=imgss[*n];
                }
	      for(i=0;i<ii;++i)
		if(wngb[i]==olab[*n])
		  break;
	      if(i==ii) 	    /* label not yet added to wngb.  */
		{
		  wngb[ii]=olab[*n];
		  ++ii;
		}
	    }
	while(++n<nf);

	/* If more than one neighboring label was found, fill in the
	   'sums' and 'nums' adjacency matrixs with the values for
	   this pixel. Recall that ii is the number of neighboring
	   labels to this river pixel. */
        if(ii>1)
          {
            rpave/=rpnum;
            for(i=0;i<ii;++i)
              for(j=0;j<ii;++j)
                if(i!=j)
                  {
                    nums[ wngb[i]* numclumps + wngb[j] ]++;
                    sums[ wngb[i]* numclumps + wngb[j] ]+=rpave;
                  }
          }
      }
  while(++ind<indf);


  /* Calculate the Signal to noise ratio of the rivers. Here, ii is
     the index in the adjacency matrix. */
  for(i=0;i<numclumps;++i)
    for(j=0;j<i;++j)
      if( nums [ ii=i*numclumps+j ] ) /* There is a connection! */
	{
	  ave = sums[ii]/nums[ii];
	  /* In case the average is nagive (only possible if sums is
	     negative), forget it. Note that even an area of 1 is
	     acceptable, and we put no area criteria here, because the
	     fact that a river exists between two clumps is
	     important. */
	  if( ave<0 )
	    {
	      nums[ii]=nums[j*numclumps+i]=0;
	      sums[ii]=sums[j*numclumps+i]=0;
	    }
	  /* Everything is ready, calculate the SN for this river:  */
	  else
	    /* Calculate the SN for this river and put it in both
	       sections of the SN adjacency matrix. Note that here we
	       want the average SN of the river, not the total. This
	       is why we haven't included the number of pixels. */
	    sns[ii] = sns[j*numclumps+i] = c * ave / sqrt(ave+err);
	}

  /* To check the matrix specific detected region:
  if(ctp->thislabel==105)
    for(i=1;i<numclumps;++i)
      {
        printf("%lu: \n", i);
        for(j=1;j<numclumps;++j)
          if(nums[ i * numclumps + j])
            printf("   %lu: %-4d %-10.2f %-10.3f\n", j,
                   nums[i*numclumps+j],
                   sums[i*numclumps+j]/nums[i*numclumps+j],
                   sns[i*numclumps+j]);
        printf("\n");
      }
  */
}





/* The true clumps were expanded. Here the job is to create two
   adjacency matrixs, which will keep the total flux on the rivers
   between two neighbouring objects. Note that since growth was not
   extended to cover the whole image, some objects might be totally
   isolated from others with no rivers existing between them. To find
   the unique neighbours of each river pixel, I will be using the same
   `wngb[]` and `ii` introduced in segmentinfo() of clumps.c.
*/
void
grownclumpstoobjects(struct clumpsthreadparams *ctp)
{
  int *nums;
  double *sums, *sns;
  long *olab=ctp->p->olab;
  size_t i, j, ii, *ind, *indf;
  size_t numclumps=ctp->numclumps;


  /* Allocate the necessary arrays: */
  errno=0; nums=calloc(numclumps*numclumps, sizeof *nums);
  if(nums==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for nums in grownclumpstoobjects "
          "(segmentation.c)", numclumps*numclumps*sizeof *nums);
  errno=0; sums=calloc(numclumps*numclumps, sizeof *sums);
  if(sums==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for sums in grownclumpstoobjects "
          "(segmentation.c)", numclumps*numclumps*sizeof *sums);
  errno=0; sns=calloc(numclumps*numclumps, sizeof *sns);
  if(sns==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for sns in grownclumpstoobjects "
          "(segmentation.c)", numclumps*numclumps*sizeof *sns);


  /* Find the signal to noise adjacency matrix for these grown
     clumps. */
  adjacencymatrixs(ctp, sns, sums, nums);


  /* Remove the connections that are smaller than the required
     threshold. Note that we don't need to worry about the number of
     pixels, because each river pixels's value was the average of its
     self and its 8 neighbors. Note that these river pixels are also
     in systematically low valued regions, which is exactly what we
     want, so they aren't completely random and any area is
     important. Connections will be based on the `nums' array. So if a
     connection is severed, its place in nums will be set to zero. */
  for(i=1;i<numclumps;++i)
    for(j=1;j<i;++j)
      if( nums[ ii=i*numclumps+j ] && sns[ii] < ctp->p->objbordersn )
	nums[ ii ] = nums[ j*numclumps+i ]=0;


  /* Find the connected regions and assign new labels (starting
     from one) to the grown clumps regions. */
  ctp->numobjects=BF_concomp_AdjMatrix(nums, numclumps, &ctp->segtoobjlabs);


  /* Correct all the pixels in olab: */
  indf = (ind=ctp->inds) + ctp->area;
  do
    if(olab[*ind]>0)
      olab[*ind]=ctp->segtoobjlabs[ olab[*ind] ];
  while(++ind<indf);


  /* Clean up */
  free(sns);
  free(nums);
  free(sums);
}



















/******************************************************************/
/*****************         Segmentation         *******************/
/******************************************************************/
void *
segmentonthread(void *inparam)
{
  struct clumpsthreadparams *ctp=(struct clumpsthreadparams *)inparam;
  struct noisechiselparams *p=ctp->p;

  float *sntable;
  char *segmentationname=p->segmentationname;
  size_t i, *ind, *indf, s0=p->lmp.s0, s1=p->lmp.s1;

  /* For the detections, there is no box, only the sides of the image. */
  ctp->x1=s0;
  ctp->y1=s1;
  ctp->x0=ctp->y0=0;

  /* Go over all the initial detections that were assigned to this
     thread. */
  for(i=0;ctp->indexs[i]!=NONTHRDINDEX;++i)
    if(ctp->indexs[i])   /* sp->indexs[i]==0 for the background. */
      {
        /* Keep the initial label of this detection. The initial label
           is mainly stored for possible debugging and checking. */
        ctp->thislabel=ctp->indexs[i];
        ctp->area=ctp->allareas[ctp->thislabel];
        ctp->inds=ctp->alllabinds[ctp->thislabel];


        /* Allocate the space for the indexs of the brightest pixels
           in each clump. We don't know how many clumps there are
           before hand, so to be safe just allocate a space equal to
           the number of pixels*/
        errno=0; ctp->topinds=malloc(ctp->area*sizeof *ctp->topinds);
        if(ctp->topinds==NULL)
          error(EXIT_FAILURE, errno, "%lu bytes for ctp->topinds in "
                "segmentonthread (segmentation.c)",
                ctp->area*sizeof *ctp->topinds);


        /* Oversegment this detection. */
        oversegment(ctp);
        if(segmentationname && p->stepnum==1)
          {free(ctp->topinds); continue;}


        /* Get the Signal to noise ratio of all the clumps within this
           detection. */
        clumpsntable(ctp, &sntable);


        /* Remove clumps with smaller signal to noise ratios and free
           all the nolonger necessary arrays: */
        removefalseclumps(ctp, sntable);
        free(sntable);
        free(ctp->xys);
        free(ctp->topinds);
        if(segmentationname && p->stepnum==2) continue;


        /* Segmenting objects can only be defined when there is more
           than one segment. Since ctp->numclumps is the number of
           clumps+1, if an object has more than one clump,
           p->numclumps will be larger than two. When there is only
           one segment (p->numseg==2) or none (p->numseg==1), then
           just set the required preliminaries to make the next steps
           be generic for all cases. */
        if(ctp->numclumps<=2)
          {
            ctp->numobjects=2;
            thisdetectionisoneobject(ctp);
            errno=0; ctp->segtoobjlabs=malloc(2*sizeof*ctp->segtoobjlabs);
            if(ctp->segtoobjlabs==NULL)
              error(EXIT_FAILURE, errno, "%lu bytes for ctp->segtoobjlabs "
                    "in segmentonthread (segmentation.c)",
                    2*sizeof*ctp->segtoobjlabs);
            if(segmentationname && ( p->stepnum>=3 || p->stepnum<=5) )
              continue;
          }
        else
          {
            /* Grow the true clumps until the growth limit. */
            prepfirstgrowth(ctp);
            growclumps(ctp, 1);
            if(segmentationname && p->stepnum==3)
              { free(ctp->blankinds); continue; }


            /* Identify the objects within the grown clumps and set
               ctp->numobjects. */
            grownclumpstoobjects(ctp);
            if(segmentationname && p->stepnum==4)
              { free(ctp->blankinds); continue; }


            /* Fill in the full detected area. */
            if(ctp->numclumps<=2)
              /* There is only one object in this whole detection. So
                 automatically, the whole region should get a label of
                 1. Do that with this function. */
              thisdetectionisoneobject(ctp);
            else
              {
                /* Grow the clumps unconditionally to fill the
                   remaining blank pixels. */
                ctp->numblanks=0; indf=(ind=ctp->inds)+ctp->area;
                do
                  if(p->olab[*ind]<0)
                    ctp->blankinds[ctp->numblanks++]=*ind;
                while(++ind<indf);
                growclumps(ctp, 0);
              }
            if(segmentationname && p->stepnum==5)
              { free(ctp->blankinds); continue; }

            /* Clean up */
            free(ctp->blankinds);
          }

        /*  */

        /* Clean up */
        free(ctp->segtoobjlabs);
      }

  /* Wait until all the threads finish: */
  if(p->cp.numthreads>1)
    pthread_barrier_wait(ctp->b);
  return NULL;
}





void
segmentdetections(struct noisechiselparams *p, size_t numobjsinit,
                  size_t *allareas, size_t **alllabinds)
{
  int err;
  pthread_t t; /* We don't use the thread id, so all are saved here. */
  pthread_attr_t attr;
  pthread_barrier_t b;
  struct clumpsthreadparams *ctp;
  size_t i, nb, *indexs, thrdcols;
  size_t numthreads=p->cp.numthreads;


  /* Allocate the arrays to keep the thread and parameters for each
     thread. */
  errno=0; ctp=malloc(numthreads*sizeof *ctp);
  if(ctp==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes in segmentdetections "
          "(segmentation.c) for ctp", numthreads*sizeof *ctp);


  /* Distribute the initial labels between all the threads.  */
  distinthreads(numobjsinit, numthreads, &indexs, &thrdcols);

  /* Spin off the threads to work on each object if more than one
     thread will be used. If not, simply start working on all the
     detections. */
  if(numthreads==1)
    {
      ctp[0].p=p;
      ctp[0].id=0;
      ctp[0].indexs=indexs;
      ctp[0].allareas=allareas;
      ctp[0].alllabinds=alllabinds;
      segmentonthread(&ctp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this running thread
	 (that spinns off the nt threads) is also a thread, so the
	 number the barrier should be one more than the number of
	 threads spinned off. */
      if(numobjsinit<numthreads) nb=numobjsinit+1;
      else                       nb=numthreads+1;
      attrbarrierinit(&attr, &b, nb);

      /* Spin off the threads: */
      for(i=0;i<numthreads;++i)
	if(indexs[i*thrdcols]!=NONTHRDINDEX)
	  {
            ctp[i].p=p;
            ctp[i].id=i;
            ctp[i].b=&b;
            ctp[i].allareas=allareas;
            ctp[i].alllabinds=alllabinds;
            ctp[i].indexs=&indexs[i*thrdcols];
	    err=pthread_create(&t, &attr, segmentonthread, &ctp[i]);
	    if(err) error(EXIT_FAILURE, 0, "Can't create thread %lu.", i);
	  }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&b);
    }

  free(ctp);
  free(indexs);
}




















/******************************************************************/
/*****************         Main function        *******************/
/******************************************************************/
void
segmentation(struct noisechiselparams *p)
{
  float *f, *fp;
  char *extname=NULL;
  long *forfits=NULL;
  size_t i, s0=p->smp.s0, s1=p->smp.s1;
  char *segmentationname=p->segmentationname;
  size_t *labareas, **labinds, numobjsinit=p->numobjects;

  /* Start off the counter for the number of objects and clumps. The
     value to these variables will be the label that is given to the
     next clump or object found. Note that we stored a copy of the
     initial number of objects in the numobjsinit variable above.*/
  p->numclumps=1;
  p->numobjects=1;


  /* Start the steps image: */
  if(segmentationname)
    {
      arraytofitsimg(segmentationname, "Input", FLOAT_IMG, p->img,
                     s0, s1, p->numblank, p->wcs, NULL, SPACK_STRING);
      arraytofitsimg(segmentationname, "Convolved-SkySubtracted",
                     FLOAT_IMG, p->conv, s0, s1, p->numblank,
                     p->wcs, NULL, SPACK_STRING);
      arraytofitsimg(segmentationname, "InitialLabels",
                     LONG_IMG, p->olab, s0, s1, 0, p->wcs,
                     NULL, SPACK_STRING);
    }

  /* All possible NaN pixels should be given the largest possible
     float flux in the convolved image (which is used for
     over-segmentation). NOTE that the convolved image is used for
     relative pixel values, not absolute ones. This is because NaN
     pixels might be in the centers of stars or bright objects (they
     might slice through a connected region). So we can't allow them
     to cut our objects. The safest approach is to start segmentation
     of each object or noise mesh with the NaNs it might contain. A
     NaN region will never be calculated in any flux measurement any
     way and if it is connecting two bright objects, they will be
     segmented because on the two sides of the NaN region, they have
     different fluxes.*/
  if(p->numblank)
    {
      fp=(f=p->conv)+s0*s1;
      do if(isnan(*f)) *f=FLT_MAX; while(++f<fp);
    }


  /* Find the true clump S/N threshold andput it in p->lmp.garray1. */
  p->b0f1=0;
  clumpsngrid(p);
  if(p->segmentationname)
    arraytofitsimg(p->segmentationname, "NoiseOversegmentaion",
                   LONG_IMG, p->clab, p->smp.s0, p->smp.s1, 0, p->wcs,
                   NULL, SPACK_STRING);



  /* Save the indexs of all the detected labels. We will be working on
     each label independently: */
  labindexs(p->olab, s0*s1, numobjsinit, &labareas, &labinds);


  /* p->clab was used for the noise clumps once, we have to reset it
     to zero. Note that olab will be initialized for each object later
     on. */
  memset(p->clab, 0, s0*s1*sizeof *p->clab);


  /* Segment the detections. When the viewer wants to check the steps,
     the operations have to be repeated multiple times to output each
     step. */
  p->b0f1=1;
  if(p->segmentationname)
    {
      p->stepnum=1;
      while(p->stepnum<6)
	{
	  memset(p->clab, 0, s0*s1*sizeof *p->clab);
          segmentdetections(p, numobjsinit, labareas, labinds);
          switch(p->stepnum)
            {
            case 1:
              extname="Over-segmentation"; forfits=p->clab; break;
            case 2:
              extname="Successful clumps"; forfits=p->clab; break;
            case 3:
              extname="Clumps grown"; forfits=p->olab; break;
            case 4:
              extname="Objects found"; forfits=p->olab; break;
            case 5:
              extname="Objects grown"; forfits=p->olab; break;
            default:
              error(EXIT_FAILURE, 0, "A bug! Please contact us at %s to "
                    "fix the problem. For some reason, the variable "
                    "p->stepnum in segmenation (segmentation.c) has the "
                    "unrecognized value of %d.", PACKAGE_BUGREPORT,
                    p->stepnum);
            }
          arraytofitsimg(p->segmentationname, extname, LONG_IMG, forfits,
                         s0, s1, 0, p->wcs, NULL, SPACK_STRING);
          ++p->stepnum;
        }
    }
  else
    segmentdetections(p, numobjsinit, labareas, labinds);


  /* Clean up */
  for(i=0;i<numobjsinit;++i) free(labinds[i]);
  free(labareas);
  free(labinds);
}