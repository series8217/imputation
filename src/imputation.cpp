//============ imputation.cpp ==========================================

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
#include <stdio.h>
#include <cmath>
#include <sys/stat.h> // for file size
#include "imputation.h"
using namespace std;

using namespace __gnu_cxx;


// Legend
// ==============================
Legend::Legend(char *fname)
{
	m_valid = readLegendFile(fname);
}

Legend::~Legend()
{
	while (m_legend.size()>0)
	{
		delete m_legend.back();
		m_legend.pop_back();
	}
}

bool Legend::readLegendFile(char *fname)
{
	
	cout << "\nOpening file " << fname << endl;
	FILE *legendFile = fopen(fname,"r");

	if (legendFile==NULL)
	{
		cout << "Failed to open legend file! Check that it exists.\n";
		return false;
	}

	cout << "Reading legend file...\n";
	// Temporary storage variables for refSNP #, first nucleotide letter, and second nucleotide letter
	rs_t rs_temp; 
	uint32_t pos_temp;
	char leg_temp[2][2]; // temp variable to store legend characters

	// Read in the legend file
	// -----------------------

	// skip the column labels
	fscanf(legendFile,"%*s%*s%*d%*d"); 

	uint32_t i;
	int read;
	i = 0;
	while (-1 != (read = fscanf(legendFile," %*2[RSrs]%d %d %[AGTC] %[AGTC]",&rs_temp,&pos_temp,leg_temp[0],leg_temp[1])))
	{
		//cout << read << " " << rs_temp << " " << leg_temp[0][0] << " " << leg_temp[1][0] << endl;
		if (read != 4)
		{
			cout << "ERROR: Malformed legend file.\n";
			return false;
		}

		// add the entry to the list
		m_legend.push_back(new LegendEntry_t);
		m_legend.back()->m_rs = rs_temp;
		m_legend.back()->m_position = pos_temp;
		m_legend.back()->m_leg[0] = leg_temp[0][0];
		m_legend.back()->m_leg[1] = leg_temp[1][0];
		m_directory[rs_temp] = m_legend.size()-1;
		i++;
	}
	fclose(legendFile);


	cout << "Finished reading legend file: " << m_legend.size() << " SNPs defined in legend." << endl;

	for (uint32_t j = 0; i < m_legend.size(); i++)
		cout << "SNP #" << i+1 << " rs" << m_legend[i]->m_rs << " 0:" << m_legend[i]->m_leg[0] << " 1:" << m_legend[i]->m_leg[1] << endl;
	
	if (i == 0)
		return false;
	
	return true;
}

SNP_t Legend::getSNP(int code, uint32_t i){
	int leg_code;
	if (code == 0 || code == 1)
	{
		SNP_t newSNP = { m_legend[i]->m_rs, m_legend[i]->m_position, m_legend[i]->m_leg[leg_code] };
		return newSNP;
	}
	else cout << "Legend::getSNP called with invalid code number " << code << endl;

	SNP_t newSNP = { 0,0,'x'};
	return newSNP;

}

SNP_t Legend::getSNP(char code, uint32_t i){
	int leg_code;
	if (code == '1') 
		leg_code = 1;
	else if (code == '0')
		leg_code = 0;
	else { cout << "Legend::getSNP called with invalid code character " << code << endl;
		SNP_t newSNP = { 0,0,'x'};
		return newSNP;
	}



	// return getSNP(leg_code, i);
	// there appears to be a kernel or compiler bug if i try to use the above line!!
	// I get a segmentation fault if I'm calling the overloaded version of this function,
	// but it goes away if i have some other instructions before it.. seems to depend on timing?
	// I don't know the specifics of the bug.

	// use these lines instead:
	SNP_t newSNP = { m_legend[i]->m_rs, m_legend[i]->m_position, m_legend[i]->m_leg[leg_code] };
	return newSNP;
	


}

uint32_t Legend::numSNPs(){
	return m_legend.size();
}



char Legend::getAllele0Char(rs_t SNP)
{
	// returns the char of the allele #0 for the given rs#
	
	return m_legend[getIndex(SNP)]->m_leg[0];
}


char Legend::getAllele1Char(rs_t SNP)
{
	// returns the char of the allele #1 for the given rs#
	
	return m_legend[getIndex(SNP)]->m_leg[1];
}

uint32_t Legend::getPosition(uint32_t i)
{
	return m_legend[i]->m_position;
}

uint32_t Legend::getPosition_rs(rs_t SNP)
{
	uint32_t index = getIndex(SNP);
	return m_legend[index]->m_position;
}

uint32_t Legend::getIndex(rs_t SNP)
{
	// returns the index in m_legend for the SNP with rs# SNP
	return m_directory[SNP];
}

rs_t Legend::getRS(uint32_t i)
{
	// returns the refSNP# for the SNP which occurs in the ith column of the input file
	return m_legend[i]->m_rs;
}

bool Legend::valid()
{
	return m_valid;
}

// returns the number code for the SNP at rs# rsnum, with the value SNPval
uint32_t Legend::getCode(rs_t rsnum, char SNPval)
{
	if (getAllele0Char(rsnum) == SNPval)
		return 0;
	else return 1;
}

// Haplotype
// ===============================


Haplotype::Haplotype(char* inphap_fname,Legend* legend)
{	// create a new haplotype from a file containing a single haplotype
	
	if (!readSingleHapFile(inphap_fname, legend))
	{
		cout << "Failed to read haplotype from single-haplotype file.\n";
		exit(1);
	}
}

Haplotype::Haplotype(FILE* hapFile,Legend* legend)
{	// create a new haplotype by reading a single line from an already-open file and legend

	if (!readSingleHap(hapFile,legend))
	{
		cout << "Failed to read haplotype from line of haplotype file.\n";
		exit(1);
	}
}



bool Haplotype::readSingleHap(FILE* hapFile, Legend* legend)
{	// reads a haplotype from a single line of an already-open file

	uint32_t i = 0; // counter of SNPs read from this haplotype (for error checking, since it should equal the # of SNPs in the legend)
	char SNPcode = 'x';
	SNP_t SNPtemp;

	int read;
	while (-1 != (read = fscanf(hapFile,"%c",&SNPcode)))
	{

		// each haplotype must be on its own line, so a newline character ends the haplotype
		if (SNPcode == '\n')
		{
			//cout << "Encountered end of haplotype.\n";
			break;
		}
		else if (SNPcode == ' ') // ignore other white space
			continue;
		else if (SNPcode == 'x') 
		{
			// 'x' indicates a missing SNP (to impute), so increment the counter but don't record a value
			i++;
			continue;
		}

		SNPtemp = legend->getSNP(SNPcode,i);
		m_SNPs[SNPtemp.m_rs] = SNPtemp.m_value;

		i++;
	}

	// if the number of SNPs read in is less than the number in the legend, the legend or input is malformed!
	if (i < legend->numSNPs())
	{
		cerr << "Error: The legend file defines " << legend->numSNPs() << " SNPs per haplotype, but the number"
		<< "of SNPs in one line of input file is " << i << "." << endl;
		return false;
	}

	return true;
}


bool Haplotype::readSingleHapFile(char* fname, Legend* legend)
{	// Reads a file with path fname containing a single haplotype
	cout << "\nOpening file " << fname << endl;
	FILE *hapFile = fopen(fname,"r");

	if (hapFile==NULL)
	{
		cout << "Failed to open input haplotype file! Check that it exists.\n";
		return false;
	}
	
	if (readSingleHap(hapFile,legend) == false)
	{
		cout << "Failed to read haplotype from file.\n";
		fclose(hapFile);
		return false;
	}
	
	fclose(hapFile);
	cout << "Finished reading input haplotype file: " << m_SNPs.size() << " tag SNPs recorded." << endl;

	return true;
}


bool Haplotype::isMissing(rs_t SNP)
{	// checks if the SNP at the given rs# has a known value in the haplotype
	char SNPval = getAllele(SNP);
	if (SNPval == 'x' || SNPval == '\0')
		return true;
	else return false;
}


char Haplotype::getAllele(rs_t SNP)
{
	// return the character of the base at SNP with rs_t SNP
	return m_SNPs[SNP];
}

uint32_t Haplotype::numSNPs()
{
	// returns the number of measured SNPs in the haplotype
	return m_SNPs.size();
}

// Correlation
// ================================
Correlation::Correlation(char *refhap_fname, Legend* legend)
{
	// makes a new correlation object. requires reference haplotype file name, and pointer to an existing legend object.
	m_legend = legend;
	m_r2threshold = 0.35;

	// set up the counts cache array
	m_counts = new uint32_t[legend->numSNPs()];
	for (uint32_t i = 0; i < legend->numSNPs(); i++)
		m_counts[i] = MAX_INT; // value of max int indicates that the count has not been made yet

	if (!readHaplotypes(refhap_fname, legend))
	{
		cout << "Failed to read haplotypes.\n";
		exit(1);
	}
}

Correlation::~Correlation()
{
	delete[] m_counts;
	
	while (m_haplotypes.size()>0)
	{
		delete m_haplotypes.back();
		m_haplotypes.pop_back();
	}

	// ********* need to iterate through the hash map of stacks and delete them all!!!!
}

bool Correlation::readHaplotypes(char *refhap_fname, Legend* legend)
{
	// reads the reference haplotype file and puts each haplotype in the m_haplotypes array
	cout << "\nOpening file " << refhap_fname << endl;
	FILE *refFile = fopen(refhap_fname,"r");

	if (refFile==NULL)
	{
		cout << "Failed to open reference haplotype file! Check that it exists.\n";
		return false;
	}
	cout << "Reading reference haplotypes file... this may take awhile.\n";

	// get the end of file position so we know when to stop reading haplotypes
	fseek(refFile,0L,SEEK_END); // go to the end of the file
	long endOfFile = ftell(refFile); // record the position
	fseek(refFile,0L,SEEK_SET); // go back to the start of the file

	while (ftell(refFile)<endOfFile) // while we haven't reached the end of the file
	{
		m_haplotypes.push_back(new Haplotype(refFile,legend)); // read a haplotype
	}

	cout << "Done reading reference haplotypes: " << m_haplotypes.size() << " haplotypes loaded from reference data.\n";

	
	fclose(refFile); // close the reference data file

	return true;

}

uint32_t Correlation::countAllele0(rs_t SNP)
{
	// counts the number of times the allele #0 of SNP1 occurs
	// on any haplotypes in the reference haplotype set

	uint32_t index = m_legend->getIndex(SNP);
	uint32_t count = m_counts[index];
	
	if (count == MAX_INT) // if the count in the array is MAX_INT, it hasn't been counted yet
	{
		char allele0char = m_legend->getAllele0Char(SNP);
	

		for (uint32_t i = 0; i< m_haplotypes.size(); i++)
		{
			if (allele0char == m_haplotypes[i]->getAllele(SNP))
				count++;
		}
		m_counts[index] = count;
	}

	return count;
}

uint32_t Correlation::countJointAllele0(rs_t SNP1, rs_t SNP2)
{
	// counts the number of times the allele #0s of SNP1 and SNP2 occur
	// on the same haplotype in the reference haplotype set
	
	uint32_t count = 0;
	char allele1char = m_legend->getAllele0Char(SNP1);
	char allele2char = m_legend->getAllele0Char(SNP2);
	

	for (uint32_t i = 0; i< m_haplotypes.size(); i++)
	{
		if ((allele1char == m_haplotypes[i]->getAllele(SNP1)) && (allele2char == m_haplotypes[i]->getAllele(SNP2)))
			count++;
	}
	
	return count;
}
    
char Correlation::getMajorAlleleChar(rs_t SNP)
{   
    // returns the character of the major allele for the SNP
    double pA0 = (double)countAllele0(SNP)/(double)numHaplotypes();

    if (pA0 >= 0.5) // if allele 0 is the major allele
	    return m_legend->getAllele0Char(SNP); // return the char for allele 0
    else return m_legend->getAllele1Char(SNP); // otherwise return the char for allele 1
}

    
char Correlation::getMinorAlleleChar(rs_t SNP)
{
    // returns the character of the minor allele for the SNP
    double pA0 = (double)countAllele0(SNP)/(double)numHaplotypes();

    if (pA0 < 0.5) // if allele 0 is the minor allele
	    return m_legend->getAllele0Char(SNP); // return the char for allele 0
    else return m_legend->getAllele1Char(SNP); // otherwise return the char for allele 1
}


int Correlation::numHaplotypes()
{
	return m_haplotypes.size();
}


float Correlation::getR(rs_t SNP1, rs_t SNP2)
{
	// returns the r correlation value for SNP1 and SNP2
	// this provides more information than the r^2 value
	// since it indicates which allele of SNP1 corresponds
	// to which allele of SNP2

	// positive r means allele 0 of SNP1 coincides with allele 1 of SNP 2
	// negative r means allele 0 of SNP2 coincides with allele 0 of SNP 2

	double pAB, pA, pB, R;
	uint32_t numHaps;
	
	numHaps = numHaplotypes();

	if (numHaps == 0)
		return 0.0;

	uint32_t countA = countAllele0(SNP1);
	uint32_t countB = countAllele0(SNP2);
	uint32_t countAB = countJointAllele0(SNP1,SNP2);

	if (countA == 0.0 || countA == 1.0)
		return 0.0;
	else pA = (double)countA/(double)numHaps;

	if (countB == 0.0 || countB == 1.0)
		return 0.0;
	else pB = (double)countB/(double)numHaps;
	

	pAB = (double)countAB/(double)numHaps;

	R = (pAB - (pA * pB)) / (sqrt(pA*(1.0-pA)) * sqrt(pB*(1.0-pB)));

	return R;
}


float Correlation::getR2(rs_t SNP1, rs_t SNP2)
{
	// returns the r^2 correlation value for SNP1 and SNP2
	double R = getR(SNP1,SNP2);
	return (R*R);
}

rs_t Correlation::getImplicator(rs_t SNP)
{
	// ***** major performance bottleneck!!! *****
	// **** need to reduce number of irrelevant implicators ***
	// **** maybe use precomputed R^2 value from hapmap?? ****

	// returns the next SNP rs# correlated with the given SNP rs#

	stack<rs_t>* impstack;
	rs_t impSNP;

	// check if a stack of implicators already exists for the given SNP
	impstack = m_implicators[SNP];

	// if m_implicators[SNP] returns a NULL pointer, the SNP is not in the table yet
	// so we need to make a new stack...
	if (impstack == NULL) 
	{
		// make a new stack of implicators
		impstack = new stack<rs_t>;
		m_implicators[SNP] = impstack; // assign it to the table
		
		// push implicating SNPs onto the stack

		// for each possible SNP
		uint32_t numSNPs = m_legend->numSNPs();
		for (uint32_t i = 0; i<numSNPs; i++)
		{
			rs_t nextSNP = m_legend->getRS(i);


			// if SNPs are more than 1,000,000 positions apart, don't use for implication
			if (m_legend->getPosition_rs(nextSNP) - m_legend->getPosition_rs(SNP) > 5000000u)
				continue;

			// if the R value is less than 0.1, skip it
			if (getR2(SNP,nextSNP) < m_r2threshold)
				continue;

			// **********************
			// should check for how many implicators there are before stopping. if less than 
			// a set amount, increase the distance threshold then find more.
			// need to try different values to see what is best for this


			impstack->push(nextSNP); // add it to the stack
		}	
		
		//cout << "\nGenerated list of " << impstack->size() << " implicators for rs" << SNP << endl;
		
	}
	
	// if any SNPs are left on the implicators stack
	if (impstack->size() > 0)
	{
		// get the top rs_t off the stack
		impSNP = impstack->top();
		impstack->pop();
	}
	else impSNP = 0; // if no SNPs left on stack return 0
			 // this should be interpreted as the indication
			 // that there are no more implicators

	return impSNP;
}

void Correlation::setR2threshold(float newthresh)
{
	if (newthresh >= 0.0 && newthresh <= 1.0)
	{
		m_r2threshold = newthresh;
	}
}
// Imputation
// ================================

// creates an imputation based on legend file, reference haplotype file, and input haplotype file
Imputation::Imputation(char *inphap_fname, char *refhap_fname, char *leg_fname)
{
	m_done = false;

	m_legend = new Legend(leg_fname);
	if (!m_legend->valid())
		exit(1);

	m_input = new Haplotype(inphap_fname,m_legend);

	m_correl = new Correlation(refhap_fname,m_legend);
}

Imputation::~Imputation()
{
	delete m_legend;

	delete m_input;
	
	delete m_correl;

	while (m_imputed.size()>0)
	{
		delete m_imputed.back();
		m_imputed.pop_back();
	}
}

void Imputation::setR2threshold(float newthresh)
{
	m_correl->setR2threshold(newthresh);
}

// runs the imputation
void Imputation::run()
{
	uint32_t numSNPs = m_legend->numSNPs();
	uint32_t numMissingSNPs = m_legend->numSNPs()-m_input->numSNPs();

	cout << "\nRunning imputation on " << numMissingSNPs << " missing SNPs... this may take awhile.\n";


	// draw a progress bar
	/*printf("[                                                  ]    0%%, 0/%d SNPs",numMissingSNPs);
    int update_every = 100;
	uint32_t bar_numBars = 50;
	float bar_Percent = 0.0;
	uint32_t bar_Progress = 0;
	char bar_bar[51];
	bar_bar[50] = 0;*/
    uint32_t imputedCount = 0;

	for (uint32_t i = 0; i < numSNPs; i++)
	{
		rs_t nextRS = m_legend->getRS(i);
		if (m_input->isMissing(nextRS))
		{
			m_imputed.push_back(impute(nextRS));
			imputedCount++;
		}
		else {
			m_imputed.push_back(new ImputedSNP_t);
			m_imputed.back()->m_SNP.m_rs = nextRS;
			m_imputed.back()->m_SNP.m_position = m_legend->getPosition_rs(nextRS);
			m_imputed.back()->m_SNP.m_value = m_input->getAllele(nextRS);
			m_imputed.back()->m_measured = true;
			m_imputed.back()->m_implicators = 0;
			m_imputed.back()->m_confidence = 1.0;
		}

        /* if (imputedCount % update_every == 0)
        {
		bar_Percent = (double)imputedCount/(double)numMissingSNPs;
		bar_Progress = bar_Percent * bar_numBars;
		uint32_t n = 0;
		for (; n<=bar_Progress; n++)
			bar_bar[n]='=';
		for (; n<bar_numBars; n++)
			bar_bar[n]=' ';
		printf("\r[%s] %.1f%%",bar_bar,bar_Percent*100);
        } */
        printf("\r%d/%d SNPs",imputedCount,numMissingSNPs);
	}

	cout << "\nFinished imputation.\n";

	m_done = true;
}

ImputedSNP_t* Imputation::impute(rs_t SNP)
{
	// imputes the value of a SNP with rs# 'SNP'

	double R;
	double R2_sum = 0.0; // accumulates positive (allele #0) and negative (allele #1) imputations (weighed R^2 values)
				                // from each implicator. the stronger the negative or positive value,
				                // the stronger the imputation in that direction.

    // Initialize the new imputed SNP object
	ImputedSNP_t* imputed = new ImputedSNP_t;
	imputed->m_SNP.m_rs = SNP;
	imputed->m_SNP.m_position = m_legend->getPosition_rs(SNP);
	imputed->m_measured = false;
	imputed->m_implicators = 0; // the number of directly measured implicating SNPs used to impute the SNP's value

	while (true) // until we break
	{
		rs_t implicator = m_correl->getImplicator(SNP); // get a correlated SNP
		if (!implicator) break; // if no correlated SNPs left, break from the while loop

		// check if a SNP was actually measured or not in the input haplotype. if it wasn't, don't impute off it.
		if (m_input->isMissing(implicator))
			continue;

		R = m_correl->getR(SNP,implicator); // get the R correlation value

		if (m_legend->getAllele0Char(implicator) == m_input->getAllele(implicator)) // if the value of the implicator in the input hap is allele 0
        {
            if (R>=0) // and the R value is positive
    			R2_sum+=R*R; // add the R^2 value
            else
                R2_sum-=R*R; // otherwise subtract the R^2 value (because this indicates allele 1 in the imputed SNP instead of allele 0
        }
		else 
        {
            if (R>=0) // if R value is positive
                R2_sum-=R*R; // subtract the R^2 value
            else 
                R2_sum+=R*R; // otherwise add the R^2 value
        }

		imputed->m_implicators++; // increment the number of implicators
	}

    // if the R^2 weighed sum is greater than or equal to zero
	if (R2_sum>0)
	{
		imputed->m_SNP.m_value = m_legend->getAllele0Char(SNP);
	}
	else if (R2_sum == 0)
    { imputed->m_SNP.m_value = m_correl->getMajorAlleleChar(SNP);
    }   
    else
    { imputed->m_SNP.m_value = m_legend->getAllele1Char(SNP);
    }

	// the "confidence" that this SNP is actually of the given value is the mean of the summed R^2 values.
    // note that this is an arbitrary measure and left to the interpretation of the user. 
    // keep in mind that a large number of low R^2 values in one direction decreases the confidence...
    // ignoring these values increases the confidence. this is NOT statistical confidence.
	if (imputed->m_implicators == 0)
		imputed->m_confidence = 0;
	else imputed->m_confidence = R2_sum / (double)imputed->m_implicators; 

	//cout << "Imputed rs" << SNP << " to " << imputed->m_SNP.m_value << " with confidence " << imputed->m_confidence << endl;

	return imputed;
}

bool Imputation::isdone()
{
	return m_done;
}

bool Imputation::printReadable(char *out_fname)
{
	if (!m_done)
	{
		cout << "print() called, but imputation not complete.\n";
		return false;
	}

	// prints the imputed haplotype to the file specified by out_fname

	cout << "\nCreating file " << out_fname << endl;

	// check if the file already exists
	FILE *outfile = fopen(out_fname,"r");
	if (outfile!=NULL)
	{
		string inBuffer;

		cout << "\n The file " << out_fname << " already exists. Overwrite? (Y/N): ";
		cin >> inBuffer;

		if (!inBuffer.compare("Y") || !inBuffer.compare("y") || !inBuffer.compare("YES")
					|| !inBuffer.compare("Yes") || !inBuffer.compare("yes"))
		{ // do nothing 
		}
		else
		{
			cout << "\nAborted printing results!\n\n";
			return false;
		}
	}

	outfile = fopen(out_fname,"w+");
	if (outfile==NULL)
	{
		cout << "\nFailed to create output haplotype file! Check that you have permission.\n";
		return false;
	}

	cout << "\nWriting imputed haplotype data...\n";
	
	// write the table header
	fprintf(outfile,"refSNP#\tposition\tSNP\tconfidence\timplicators\n");

	for (uint32_t i = 0; i < m_imputed.size(); i++)
	{
		fprintf(outfile,"rs%d\t%d\t%c\t%.3f\t%d\n",m_imputed[i]->m_SNP.m_rs, m_imputed[i]->m_SNP.m_position,
							m_imputed[i]->m_SNP.m_value, m_imputed[i]->m_confidence,
							m_imputed[i]->m_implicators);
	}

	return true;
}


bool Imputation::printRaw(char *out_fname)
{
	if (!m_done)
	{
		cout << "printRaw() called, but imputation not complete.\n";
		return false;
	}

	// prints the imputed haplotype to the file specified by out_fname, in raw format

	cout << "\nCreating file " << out_fname << endl;

	// check if the file already exists
	FILE *outfile = fopen(out_fname,"r");
	if (outfile!=NULL)
	{
		string inBuffer;

		cout << "\n The file " << out_fname << " already exists. Overwrite? (Y/N): ";
		cin >> inBuffer;

		if (!inBuffer.compare("Y") || !inBuffer.compare("y") || !inBuffer.compare("YES")
					|| !inBuffer.compare("Yes") || !inBuffer.compare("yes"))
		{ // do nothing 
		}
		else
		{
			cout << "\nAborted printing results!\n\n";
			return false;
		}
	}

	outfile = fopen(out_fname,"w+");
	if (outfile==NULL)
	{
		cout << "\nFailed to create output haplotype file! Check that you have permission.\n";
		return false;
	}

	cout << "\nWriting raw imputed haplotype data...\n";
	
	for (uint32_t i = 0; i < m_imputed.size(); i++)
	{
		fprintf(outfile,"%d ", m_legend->getCode(m_imputed[i]->m_SNP.m_rs,m_imputed[i]->m_SNP.m_value));
	}
	fprintf(outfile," \n");

	fclose(outfile);
	cout << "Finished writing imputed haplotype data.\n";

	return true;
}
