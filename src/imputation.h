#ifndef __IMPUTATION_INCLUDED___
#define __IMPUTATION_INCLUDED___
//============ imputation.h ==========================================


// Imputation software for CS124 Computational Genetics, Spring 2008, UCLA

// by Steven Thomas Snyder, stsnyder@ucla.edu




/*(C) copyright 2008, Steven Snyder, All Rights Reserved

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
*/


#include <string>
#include <ext/hash_map>
#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <stdint.h>
using namespace std;

using namespace __gnu_cxx; // for ext/hash_map

const uint32_t MAX_INT=4294967295u; // 2^32 -1

// SNP_t
// =============================
typedef uint32_t rs_t; // the type for the reference value/name/whatever of a SNP

struct eqrs
{
	bool operator()(const rs_t rs1, const rs_t rs2) const
	{
		return (rs1 == rs2);
	}
};

// SNP_t
// ============================

// Used for passing around complete SNPs with their refSNP #.
struct SNP_t
{

	rs_t m_rs;    // RefSNP
	uint32_t m_position; // position
	char m_value; // A, G, T, C, x, the nucleotide at this location (or x if unknown)
};


// ImputedSNP_t
// ========================

// contains information about a missing SNP filled in by imputation
struct ImputedSNP_t
{
	SNP_t m_SNP;
	bool m_measured; // true if the SNP was actually directly measured. in this case the 
			 // members of ImputedSNP_t below are arbitrary and/or undefined.
			// this is used because in some cases we want to store directly measured
			// SNPs with imputed SNPs, such as in the final haplotype of mixed
			// imputations and measured SNPs.

	float m_confidence; // the confidence that this SNP is actually of the given value.
				// this is calculated by multiplying the power of the correlation/association
				// test by the probability of the SNP having the given value assuming
				// the correlation is correct.

	uint32_t m_implicators; // the number of directly measured implicating SNPs used to impute the SNP's value
	
};


// Legend
// ==============================
class Legend
{
private:
	/* The Legend class contains all the information and methods necessary to parse
	 * a legend file. */
	
	struct LegendEntry_t {
		/* A LegendEntry_t must exist for each SNP in the haplotypes input file.
		 * Each entry stores the reference number and character translation for
		 * the haplotypes input file. */
		rs_t m_rs; // RefSNP number
		char m_leg[2];	// m_leg[0] is the nucleotide letter for '0' in the haplotypes file.
				// m_leg[1] is the nucleotide letter for '1' in the haplotypes file.
		uint32_t m_position; // SNP position
	};
	
	vector<LegendEntry_t*> m_legend; // pointers to the legend entires

	hash_map <rs_t,uint32_t> m_directory; // lookup map of refSNP rs# to legend entry index in m_legend

	bool m_valid; // true if the legend was formed properly at construction, false otherwise

public:
	// Creates a new legend from a filename for the legend file
	Legend(char *fname);

	~Legend();

	// Reads a legend file in as the legend
	bool readLegendFile(char *fname);

	// returns the SNP for code char ('0' or '1')  where the SNP occurs in the ith column of the input file
	SNP_t getSNP(char code, uint32_t i);

	// returns the SNP for code integer (0 or 1), where the SNP occurs in the ith column of the input file
	SNP_t getSNP(int code, uint32_t i);

	// returns the refSNP# for the SNP which occurs in the ith column of the input file
	rs_t getRS(uint32_t i);
	
	// returns the position of a SNP by index #
	uint32_t getPosition(uint32_t i);

	// returns the position of a SNP by rs #
	uint32_t getPosition_rs(rs_t SNP);

	// returns the char of the allele #0 for the given rs#
	char getAllele0Char(rs_t SNP);

	// returns the char of the allele #1 for the given rs#
	char getAllele1Char(rs_t SNP);

	// returns the index in m_legend for the SNP with rs# SNP
	uint32_t getIndex(rs_t SNP);

	// returns the number code for the SNP at rs# rsnum, with the value SNPval
	uint32_t getCode(rs_t rsnum, char SNPval);

	// returns the # of SNPs in the legend
	uint32_t numSNPs();

	bool valid();
};


// Haplotype
// ===============================
class Haplotype
{
private:
	/* The Haplotdelete arrayype class contains methods for parsing and storing haplotype data,
	 * given a Legend for the haplotype file(s) */
	hash_map <rs_t,char> m_SNPs; // the haplotype: base pair letters in a hash_map indexed by refSNP rs#
					// note that the value of the hash_map is not SNP_t, it is a character!!

public:
	// create a new haplotype from a file containing a single haplotype, and an already-open file
	Haplotype(char* inphap_fname,Legend* legend);

	// create a new haplotype by reading a single line from an already-open file and legend
	Haplotype(FILE* hapFile,Legend* legend);

	// returns the number of SNPs in the haplotype
	uint32_t numSNPs();

	// read one line from a file containing possibly many haplotypes
	bool readSingleHap(FILE* hapFile, Legend* legend);

	// checks if the SNP at the given rs# has a known value in the haplotype
	bool isMissing(rs_t SNP);

	// read a file containing a single haplotype
	bool readSingleHapFile(char* inphap_fname,Legend* legend);

	// return the character of the base at SNP with rs_t SNP
	char getAllele(rs_t SNP);

};

// Correlation
// ================================
class Correlation
{
private:
	vector<Haplotype*> m_haplotypes; // storage for reference haplotypes
	Legend* m_legend; // pointer to legend for reference haplotypes

	hash_map <rs_t, stack<rs_t>* > m_implicators; // map to stacks for implicators of SNPs.

	uint32_t *m_counts; // cache of counts of alleles

	double m_r2threshold; // minimum r^2 value for a SNP to be used as an implicator in imputation

	// *CAN BE IMPROVED*
	// Right now all correlations are generated on the fly. N imputations will take N*t
	// where t is the time required to do a single imputation. By caching the correlations
	// this would go a lot faster!
	// More importantly, missing SNPs correlated with other missing SNPs cannot be imputed
	// without multiple passes through the haplotype and this may not be very accurate...

	// need to figure out a data type to use for storing the correlations. 
	// should keep track of distance between SNPs so we dont cache correlations for far-away
	// SNPs as we move down the haplotype during imputation..
	// a complete 4-byte-per-datum correlation matrix would be 150 gigagbytes!

public:

	Correlation(char *refhap_fname, Legend* legend);

	~Correlation();

	// reads the reference haplotypes file
	bool readHaplotypes(char * refhap_fname, Legend* legend);


	// returns the r correlation value for SNP1 and SNP2
	// this provides more information than the r^2 value
	// since it indicates which allele of SNP1 corresponds
	// to which allele of SNP2
	float getR(rs_t SNP1, rs_t SNP2);

	// returns the r^2 correlation value for SNP1 and SNP2
	float getR2(rs_t SNP1, rs_t SNP2);

	// returns the number of haplotypes in the correlation set
	int numHaplotypes();

	// counts the number of times the allele #0 of SNP1 occurs
	// on any haplotypes in the reference haplotype set
	uint32_t countAllele0(rs_t SNP);

	// counts the number of times the alleles #0 of SNP1 and SNP2 occur
	// on the same haplotype in the reference haplotype set
	uint32_t countJointAllele0(rs_t SNP1, rs_t SNP2);

	// returns the next SNP that the given SNP is correlated with
	rs_t getImplicator(rs_t SNP);

    // returns the character of the major allele for the SNP
    char getMajorAlleleChar(rs_t SNP);

    // returns the character of the minor allele for the SNP
    char getMinorAlleleChar(rs_t SNP);

	void setR2threshold(float newthresh);


};

// Imputation
// ================================
// This is the "master" class for running an imputation. Call the constructor
// and then the run() method on the newly generated Imputation. isdone() will return
// true when the imputation is done running. 
// To get the data, use printImputation() with the appropriate arguments
// (file descriptor, etc) for where you want the imputation data to be written to.
class Imputation
{
private:
	Legend* m_legend;
	Haplotype* m_input;
	Correlation* m_correl;

	vector<ImputedSNP_t*> m_imputed; // ***** don't forget to free these!!!

	bool m_done;

public:
	// creates an imputation based on legend file, reference haplotype file, and input haplotype file
	Imputation(char *inphap_fname, char *refhap_fname, char *leg_fname);

	~Imputation();

	// runs the imputation
	void run();

	// returns true if the imputation has been run and finished successfully
	bool isdone();

	// imputes the value of SNP
	ImputedSNP_t* impute(rs_t SNP);

	// prints the results of the imputation (if done) in readable form
	bool printReadable(char *out_fname);

	// prints the results of the imputation (if done) in the same format as the input file
	bool printRaw(char *out_fname);

	// set the minimum r^2 value for a SNP to be used as an implicator in imputation
	void setR2threshold(float newthresh);

};

#endif //include guard
