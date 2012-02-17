(C) copyright 2008, Steven Snyder, All Rights Reserved

LICENSING INFORMATION:
 This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Imputation software by Steven Snyder <stsnyder@ucla.edu>
June 4, 2008
Version 0.62


Use:  imputation INPUT REFERENCE LEGEND RAWOUTFILE [DETAILOUTFILE]

Where:
	INPUT is the input haplotype to be imputed on

	REFERENCE is the input haplotypes to use for association statistics

	LEGEND is the file specifying the rs#, position, and nucleotide for
	       the '0's and '1's in the incomplete data and haplotypes input files.
	       these should all be in the same format as the HapMap phased
	       haplotype data.

	RAWOUTFILE is the file the imputed haplotype should be written to
	           in the same format as the input file

	[DETAILOUTFILE] (optional) is the file the imputed haplotype should be
	          written to in detailed, readable form



Notes on input files:

INPUT - This file must have a single haplotype in the same format as the
		   reference haplotypes file (see below)

REFERENCE - This file should be in the same format as the HapMap phased
		haplotype files. Example for 16-SNP haplotypes:
		
		0 0 1 1 0 1 0 1 0 0 0 0 1 1 0 1
		1 1 0 0 1 1 1 0 1 0 1 0 1 1 0 1
		1 1 1 1 0 0 1 1 1 1 1 1 0 0 1 1
		1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0
		0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1


		There may be any number of haplotypes. Make sure that the last character is
		the last base number code for the last haplotype. Do not put a newline/return
		after the last haplotype!

LEGEND - The legend file must have the following format:

	rs	position	0	1
	rs[rs#]	[pos#]	    [A,G,T,C]  [A,G,T,C]
	rs[rs#]	[pos#]	    [A,G,T,C]  [A,G,T,C]
	rs[rs#]	[pos#]	    [A,G,T,C]  [A,G,T,C]
	rs[rs#]	[pos#]	    [A,G,T,C]  [A,G,T,C]
	 ...      ...          ...       ...
	rs[rs#]	[pos#]	    [A,G,T,C]  [A,G,T,C]

	Note that the rs#s should be preceded by "rs" (without the quotes) as in the HapMap
	phased haplotype files available from hapmap.org.



