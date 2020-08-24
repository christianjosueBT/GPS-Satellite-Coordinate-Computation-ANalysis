// Summary:
//    Contains the NL_NRinexUtils_H namespace and several functions for working 
//    with RINEX data.
//
// History:
//    Jan 09/13 - Created by Mark Petovello
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

#ifndef NL_NRinexUtils_H
#define NL_NRinexUtils_H

#include <string>
#include <iostream>

#include "rinex.h"

namespace NRinexUtils
{
   //**
   // Summary:
   //    Get the type of RINEX file (e.g., observations, navigation, atc.)
   //
   // Arguments:
   //    None.
   //
   // Returns:
   //    A char containing the type of file.  This can be one of the typical characters that 
   //    end the file name (e.g., 'O', 'N', 'G', 'M', etc.)
   char GetFileType( const std::string& rinexFilename );


   //**
   // Summary:
   //    Open a generic RINEX file for input.
   //
   // Arguments:
   //    rinexFile - The file to be opened.
   //    rinexFilename - The name of the file.
   //
   // Returns:
   //    True if successful and false otherwise.
   bool OpenRinexFileForInput( NGSrinex::RinexFile &rinexFile, const std::string& rinexFilename );


   //**
   // Summary:
   //    Open a RINEX observation file for input.
   //
   // Arguments:
   //    rinexFile - The file to be opened.
   //    rinexFilename - The name of the file.
   //
   // Returns:
   //    True if successful and false otherwise.
   bool OpenRinexObservationFileForInput( NGSrinex::RinexObsFile &rinexFile, const std::string& rinexFilename );


   //**
   // Summary:
   //    Open a RINEX navigation (GPS ephemeris) file for input.
   //
   // Arguments:
   //    rinexFile - The file to be opened.
   //    rinexFilename - The name of the file.
   //
   // Returns:
   //    True if successful and false otherwise.
   bool OpenRinexNavigationFileForInput( NGSrinex::RinexNavFile &rinexFile, const std::string& rinexFilename );

   void RinexGpsEphemeris( const std::string& rinexFilename, std::ostream& binaryOutFile );

   //**
   // Summary:
   //    Open a RINEX observation file for output.
   //
   // Arguments:
   //    rinexFile - The file to be opened.
   //    rinexFilename - The name of the file.
   //
   // Returns:
   //    True if successful and false otherwise.
   bool OpenRinexFileForOutput( NGSrinex::RinexFile &rinexFile, const std::string& rinexFilename );
   
   struct stGpsEphemeris
      {
         //**
         // Satellite group delay in units of seconds.
         double         tgd;

         //**
         // Clock issue data (unitless).
         short int         aodc;

         //**
         // Reference time of clock in units of seconds into the week.
         double         toc;

         //**
         // Satellite clock drift rate parameter in units of second per second-squared.
         double         af2;

         //**
         // Satellite clock drift parameter in units of seconds per seconds.
         double         af1;

         //**
         // Satellite clock offset parameter in units of seconds.
         double         af0;

         //**
         // Orbit issue data (unitless).
         unsigned char         aode;

         //**
         // Amplitude of sine correction to orbital radius in units of metres.
         float        crs;

         //**
         // Mean anomaly correction in units of radians per second.
         double         deltan;

         //**
         // Mean anomaly at reference time in units of radians.
         double         m0;

         //**
         // Amplitude of cosine correction to argument of latitutde in units of radians.
         double         cuc;

         //**
         // Orbital eccentricity (unitless).
         double         e;

         //**
         // Amplitude of sine correction to argument of latitutde in units of radians.
         double         cus;

         //**
         // Square root of semi-major axis of orbit in units of sqrt(metres).
         double         sqrta;

         //**
         // Reference time of orbit in units of seconds into the week.
         double        toe;

         //**
         // Amplitude of cosine correction to inclination angle in units of radians.
         double         cic;

         //**
         // Longitude of ascending node at reference time in units of radians.
         double         omega0;

         //**
         // Amplitude of sine correction to inclination angle in units of radians.
         double         cis;

         //**
         // Inclination angle at reference time in units of radians.
         double         i0;

         //**
         // Amplitude of cosine correction to orbital radius in units of metres.
         float         crc;

         //**
         // Argument of perigee in units of radians.
         double         w;

         //**
         // Rate of right ascension in units of radians per second.
         double         omegadot;

         //**
         // Rate of inclination in units of radians per second.
         double         idot;

      }; // end struct stGpsEphemeris

};


#endif //NL_NRinexUtils_H