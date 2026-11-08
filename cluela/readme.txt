cluela - Data Clustering by k-means with extra features:
     supports numeric and category data
     supports missing data
     can accept text labels as input
     k-medians available
     can measure agreement with a known partition
     algorithms faster than O(NPK)

address correspondence to:  andy@13olive.net
See the comments in the source code for usage instructions.
An article discussing the package is available at www.codeproject.com

************************************************************************
                              CONTENTS
************************************************************************
eigen.f       Solve symmetric eigenproblems, by Martin, et. al.
eigen.c       Public domain

sort.f        Sort lists and apply permutations, by R. Singleton, et al.
sort.c        Public domain

kmns.f        Data clustering by the k-means algorithm
kmns.c        by J. A. Hartigan and M. A. Wong, Yale University.
              Copyright (c) 1979, Royal Statistical Society
              This program may be used and copied provided no fee is charged.

median.f      The median of an array, by Alan Miller, CSIRO
median.c      Public domain

safdiv.f      Test if division is safe, based on a program by Julio.
safdiv.c      Public domain

fnv32.f       The FNV-1a hash function, by Landon Curt Noll
fnv32.c       Public domain

stats.f       Functions for data clustering, by Andy Allinger
stats.c       Public domain
prep.f
prep.c
osiman.f
osiman.c
mds.f
mds.c
ght.f
ght.c
cluela.f
cluela.c
cluste.f
cluste.c
neugas.f
neugas.c
valid.f
valid.c
dump.f
dump.c

machine.f     Utility functions in the public domain.
random.f
array.c

t_crx.f       Test programs in the public domain.
t_crx.c
t_embed.f     
t_embed.c
t_fire.f
t_fire.c
t_horse.f
t_horse.c
t_round.f
t_round.c
t_soy.f
t_soy.c
t_votes.f
t_votes.c
t_wine.f
t_wine.c

crx.data                          Data files by various authors.
forestfires.csv                   Distributed FOR RESEARCH PURPOSES ONLY!
horse-colic.both 
house-votes-84.data
soybean.both    
wine.data

article.pdf              A discussion of data clustering in general
                         and this package in particular
readme.txt               This file

************************************************************************
                        IMPORTANT DISCLAIMER
************************************************************************
  The software is intended for research and educational purposes.  While
no restrictions are placed on its use, there is ABSOLUTELY NO WARRANTY,
neither by the authors nor by codeproject.com
  The names and affiliations of authors of subroutines included in the
package are for acknowledgment only, and in no way indicates their
endorsement of this software distribution.
  This software MUST NOT be used in applications where life or property
is at stake until the user has investigated whether it meets his needs.
************************************************************************

************************************************************************
                         Revision History
************************************************************************
3   24 Dec 2020   Calling sequences have changed 
                  Delete invalid missing-data treatment from PCA
                  Add multi-dimensional scaling (MDS)
                  Maximum distance criterion in binary tree search
                  Removed arithmetic IF, EQUIVALENCE and union from binary tree
                  Replacement for safe-to-divide
                  Remove compensated addition
                  Use intrinsic string comparison in IREID
                  Standardize form of distance function to weighted average
                  Sorting in CENT only when robust is true.
                  Forest fire problem replaces Maple Leaf Rag problem.
                  New error #11 for weights that sum to zero

2    8 Nov 2017   Fixed mistake in USORTU that returned a bad index
                  when N equals 1.

1    7 Oct 2017   Changed loop in KCLUST for deleting clusters to
                  negative step.  This should fix mistakes that could
                  occur when more than one cluster is to be deleted.

0   27 Aug 2017
************************* end of readme.txt ****************************
