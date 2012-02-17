//============ main.cpp ==========================================


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
#include "imputation.h"

using namespace std;

using namespace __gnu_cxx;


// function void printHelp() ---- prints the command line option information to the console
// ==================================
void printHelp(){
	cout << " Use:\n\n  imputation INPUT REFERENCE LEGEND RAWOUTFILE [DETAILOUTFILE]\n\n";
	cout << " Where:\n\n";
	cout << "  INPUT is the input haplotype to be imputed on\n\n";
	cout << "  REFERENCE is the input haplotypes to use for association statistics\n\n";
	cout << "  LEGEND is the file specifying the rs#, position, and nucleotide for\n";
	cout << "          the '0's and '1's in the incomplete data and haplotypes input files.\n";
	cout << "          these should all be in the same format as the HapMap phased\n";
	cout << "          haplotype data.\n\n";
	cout << "  RAWOUTFILE is the file the imputed haplotype should be written to\n";
	cout << "          in the same format as the input file\n\n";
	cout << "  [DETAILOUTFILE] (optional) is the file the imputed haplotype should be\n";
	cout << "          written to in detailed, readable form\n\n";
}


// function main
// ===================================
int main (int argc, char **argv)
{
	Imputation *imputation;
	char *incompDataFilename;
	char *haplotypesFilename;
	char *legendFilename;
	char *outRawFilename;
	char *outDetailFilename;
	float r2threshold;
	string inBuffer;

	cout << "\n";
	
	if (argc != 5 && argc !=6 && argc !=7)
	{	cout << "Invalid command line options.\n\n";
		printHelp();
		exit(1);
	}
	
	incompDataFilename = argv[1];

	haplotypesFilename = argv[2];

	legendFilename = argv[3];

	outRawFilename = argv[4];

	if (argc >= 6)
		outDetailFilename = argv[5];
	else outDetailFilename = NULL;

	if (argc == 7)
		sscanf(argv[6],"%f",&r2threshold);
	else r2threshold = 0.35;

	

	cout << "Using the following files for input data:\n\n";
	cout << " Incomplete data set: " << incompDataFilename << endl;
	cout << "Reference haplotypes: " << haplotypesFilename << endl;
	cout << "              Legend: " << legendFilename << endl;
	cout << "     Raw Output File: " << outRawFilename << endl;
	cout << "Readable Output File: ";
	if (outDetailFilename)
		cout << outDetailFilename << endl;
	else cout << "N/A" << endl;
	cout << "	r^2 threshold: " << r2threshold << endl;


	cout << "\n Continue? (Y/N): ";
	cin >> inBuffer;
	//cout << "Y\n";
	//inBuffer = "Y";

	if (!inBuffer.compare("Y") || !inBuffer.compare("y") || !inBuffer.compare("YES")
				|| !inBuffer.compare("Yes") || !inBuffer.compare("yes"))
	{
		// run a new imputation
		imputation = new Imputation(incompDataFilename,haplotypesFilename,legendFilename);
		imputation->setR2threshold(r2threshold);
		imputation->run();

		imputation->printRaw(outRawFilename);

		if (outDetailFilename)
			imputation->printReadable(outDetailFilename);

		if (imputation->isdone())
			cout << "\nFinished imputation.\n\n";

        cout << "\n Clearing temporary files.\n";

		delete imputation;

        cout << "\n Done! \n";

	}
	else cout << "\nAborted!\n\n";



	return 0;
}
