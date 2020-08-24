// Summary:
//    Skeleton program to open and read a RINEX ephemeris (navigation) file.  Students in 
//    ENGO 465 will use this to develop code to compute satellite positions at specified times.
//
// History:
//    Feb 04/13 - Created by Mark Petovello
//
// Copyright:
//    Position, Location And Navigation (PLAN) Group
//    Department of Geomatics Engineering
//    Schulich School of Engineering
//    University of Calgary
//
// Disclaimer:
//    This source code is not freeware nor shareware and is only provided under 
//    an agreement between authorized users/licensees and the University of 
//    Calgary (Position, Location And Navigation (PLAN) Group, Geomatics 
//    Engineering, Schulich School of Engineering) and may only be used under 
//    the terms and conditions set forth therein.

#include <iostream>
#include <fstream>
#include <iomanip>

#include "rinex.h"
#include "NRinexUtils.h"
#include <Eigen/Dense>
#include <vector>
#include <cmath>

using namespace std;
using namespace NGSrinex;
using namespace Eigen;

#define M_PI 3.1415926535898
int main(int argc, char* argv[]) {

	// input file names
	string rinexFilename;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// START: Get ephemeris data

	rinexFilename = "brdc0070.20n";

	// open the input RINEX file
	RinexNavFile rinexInput;
	if (!NRinexUtils::OpenRinexNavigationFileForInput(rinexInput, rinexFilename))
	{
		cout << "NRinexUtils::OpenRinexNavigationFileForInput() - Could not open RINEX file." << endl;
		return 0;
	}

	// END: Get ephemeris data
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// START: Open Output Files

	ofstream outFile;
	outFile.open("brdcbest.txt");
	outFile << fixed << setprecision(5);
	ofstream outFile2;
	outFile2.open("brdc24hr.txt");
	outFile2 << fixed << setprecision(5);

	if (outFile.fail())
	{
		cout << "Could not open brdcbest.xyz file..." << endl;
		exit(1);
	}
	if (outFile2.fail())
	{
		cout << "Could not open brdc24hr.xyz file..." << endl;
		exit(1);
	}

	// End: Open Output Files
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// START: Read the ephemeris data

	const double mu = 3.986005E14;
	const double omegaE = 7.292115147E-5;

	// GPS ephemeris data from RINEX file
	NGSrinex::PRNBlock  currentRinexEph;
	vector<NGSrinex::PRNBlock> PRN11;    

   // read the ephemeris
   try
   {
      // for each set of ephemeris parameters in the RINEX file...
      while( rinexInput.readPRNBlock( currentRinexEph ) != 0 )
      {
         // check for unhealthy satellite
		  if (currentRinexEph.getSvHealth() != 0)
			  continue;

		   double PRN = currentRinexEph.getSatellitePRN();

		   if (PRN != 11)
			   continue;

		   PRN11.push_back(currentRinexEph);

	  }//for each set of ephemeris parameters in the RINEX file...

	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // START: TASK 2.1 and 2.2

	  double t = 0;

	  for (int i = 0; i < 96; i++) {

		  int location = 0;
		  currentRinexEph = PRN11.at(0);
		  t = 172800 + i * 900;
		  double condition = abs(t - currentRinexEph.getToe());

		  for (int j = 0; j < PRN11.size(); j++) {
			  double x = abs(t - PRN11.at(j).getToe());
			  if (x < condition) {
				  location = j;
				  condition = x;
			  }
		  }
		  //Set the current RINEX parameters to be the closest set for all positions
		  currentRinexEph = PRN11.at(location);
		  double toe = currentRinexEph.getToe();
		  double tk = t - toe;

		  // Step 1: Compute mean anomaly and mean motion
		  double a = pow(currentRinexEph.getSqrtA(), 2);    
		  double n0 = sqrt(mu / pow(a, 3));
		  double n = n0 + currentRinexEph.getDeltan();
		  double Mk = currentRinexEph.getMo() + n * tk;  

		  // Ensure Mk is between 0 and 2Pi
		  while (Mk >= 2.0 * M_PI || Mk < 0.0) {
			  if (Mk >= 2.0 * M_PI)
				  Mk -= 2.0 * M_PI;
			  if (Mk < 0)
				  Mk += 2.0 * M_PI;
		  }

		  // Step 2: Solve eccentric anomaly iteratively
		  double Ek = Mk;
		  for (int k = 0; k < 6; k++) {
			  Ek = Mk + currentRinexEph.getEccen() * sin(Ek);
		  }

		  // Step 3: Compute true anomaly
		  double vK = atan2(sqrt(1 - pow(currentRinexEph.getEccen(), 2)) * sin(Ek), cos(Ek) - currentRinexEph.getEccen(););

		  // Step 4: Compute argument of latitude of the satellite
		  double latK = currentRinexEph.getLilOmega() + vK;

			// Step 5: Compute corrections to Keplerian orbit
		  double delta_uK = currentRinexEph.getCuc() * cos(2 * latK) +
			  currentRinexEph.getCus() * sin(2 * latK);
		  double delta_rK = currentRinexEph.getCrc() * cos(2 * latK) +
			  currentRinexEph.getCrs() * sin(2 * latK);
		  double delta_iK = currentRinexEph.getCic() * cos(2 * latK) +
			  currentRinexEph.getCis() * sin(2 * latK);

			// Step 6: Compute corrected values
		  double uK = latK + delta_uK;
		  double rK = a * (1 - currentRinexEph.getEccen() * cos(Ek)) + delta_rK;
		  double iK = currentRinexEph.getIo() + currentRinexEph.getIdot() * tk + delta_iK;

			// Step 7: Correct longitude of ascending node
		  double omegaK = (currentRinexEph.getBigOmega() - omegaE * toe) +
			  (currentRinexEph.getBigOmegaDot() - omegaE) * tk;
		  omegaK = fmod(omegaK, 2 * M_PI);

			 // Step 8: Compute position of the GPS satellite in the orbital plane
		  double x = rK * cos(uK);
		  double y = rK * sin(uK);

		  // Step 9: Calculate satellite coordinate in e frame
		  double xSV = x * cos(omegaK) - y * sin(omegaK) * cos(iK);
		  double ySV = x * sin(omegaK) + y * cos(omegaK) * cos(iK);
		  double zSV = y * sin(iK);

		  //Output the values to the 24-hour file, "brdcbest.xyz"
		  outFile << t << "   " << xSV << "   " << ySV << "   " << zSV << endl;

		  //Update time
		  t = t + 900;
	  }

	  // ~~~~~~~~
	  // TASK 2.3
	  // ~~~~~~~~

	  t = 0;

	  // Find the 24-hour ephemeris
	  for (int i = 0; i < 96; i++) {

		  //Setting the current RINEX parameters to be the first set for all positions
		  currentRinexEph = PRN11.at(0);
		  double Toe = currentRinexEph.getToe();
		  t = 172800 + i * 900;

		  // Step 1: Compute mean anomaly and mean motion
		  double tK = t - Toe;
		  double A = pow(currentRinexEph.getSqrtA(), 2);
		  double n0 = sqrt(mu / pow(A, 3));
		  double n = n0 + currentRinexEph.getDeltan(); 
		  double Mk = currentRinexEph.getMo() + n * tK;										 

		  // Ensure Mk is between 0 and 2Pi
		  while (Mk >= 2.0 * M_PI || Mk < 0.0) {
			  if (Mk >= 2.0 * M_PI)
				  Mk -= 2.0 * M_PI;
			  if (Mk < 0)
				  Mk += 2.0 * M_PI;
		  }

		  // Step 2: Solve eccentric anomaly iteratively
		  double Ek = Mk;
		  for (int k = 0; k < 6; k++) {
			  Ek = Mk + currentRinexEph.getEccen() * sin(Ek);
		  }

		  // Step 3: Compute true anomaly
		  double vK = atan2(sqrt(1 - pow(currentRinexEph.getEccen(), 2)) * sin(Ek), cos(Ek) - currentRinexEph.getEccen());

		  // Step 4: Compute argument of latitude of the satellite
		  double latK = currentRinexEph.getLilOmega() + vK;

			// Step 5: Compute corrections to Keplerian orbit
		  double delta_uK = currentRinexEph.getCuc() * cos(2 * latK) +
			  currentRinexEph.getCus() * sin(2 * latK); 
		  double delta_rK = currentRinexEph.getCrc() * cos(2 * latK) +
			  currentRinexEph.getCrs() * sin(2 * latK);
		  double delta_iK = currentRinexEph.getCic() * cos(2 * latK) +
			  currentRinexEph.getCis() * sin(2 * latK);

			// Step 6: Compute corrected values
		  double uK = latK + delta_uK;
		  double rK = A * (1 - currentRinexEph.getEccen() * cos(Ek)) + delta_rK;
		  double iK = currentRinexEph.getIo() + currentRinexEph.getIdot() * tK + delta_iK;

			// Step 7: Correct longitude of ascending node
		  double omegaK = (currentRinexEph.getBigOmega() - omegaE * Toe) +
			  (currentRinexEph.getBigOmegaDot() - omegaE) * tK;
		  omegaK = fmod(omegaK, 2 * M_PI);

		  // Step 8: Compute position of the GPS satellite in the orbital plane
		  double x = rK * cos(uK);
		  double y = rK * sin(uK);

		  // Step 9: Calculate satellite coordinate in e frame
		  double xSV = x * cos(omegaK) - y * sin(omegaK) * cos(iK);
		  double ySV = x * sin(omegaK) + y * cos(omegaK) * cos(iK);
		  double zSV = y * sin(iK);

		  //Output the values to the 24-hour file, "brdc24hr.xyz"
		  outFile2 << t << "   " << xSV << "   " << ySV << "   " << zSV << endl;

		  //Update time
		  t = t + 900;
	  }

   }
   catch (RinexReadingException & readingExcep)
   {
	   cout << " RinexReadingException is: " << endl << readingExcep.getMessage() << endl;
   }
   // END: Read the ephemeris data

   // close files (disconnect the stream from the file)
   outFile.close();
   outFile2.close();

   // end of program
   system("pause");

   return 0;
}
