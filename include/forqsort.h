/*********************************************************************
forqsort -- Functions used by qsort to sort an array.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef FORQSORT_H
#define FORQSORT_H

/* Pointer used to sort the indexs of an array based on their flux
   (value in this array). */
extern float *forqsortindexarr;

int
indexfloatdecreasing(const void * a, const void * b);

int
intdecreasing(const void * a, const void * b);

int
intincreasing(const void * a, const void * b);

int
floatdecreasing(const void * a, const void * b);

int
floatincreasing(const void * a, const void * b);

int
doubledecreasing(const void * a, const void * b);

int
doubleincreasing(const void * a, const void * b);

#endif