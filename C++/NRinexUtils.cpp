// Summary:
//    Contains the implementation of functions in the NL_NRinexUtils_H namespace.
//
// History:
//    Jan 09/13 - Created by Mark Petovello
//    Jan 24/14 - Mark Petovello updated GetFileType() to use a default constructor instead of a 
//                copy constructor when creating a RinexFile object.  This is necessary for compatibility
//                with compilers other than Microsoft Visual Studio.
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

#include "NRinexUtils.h"

using namespace NGSrinex;

char NRinexUtils::GetFileType( const std::string& rinexFilename )
{
   // create a RINEX file object
   RinexFile rinexFile;

   try 
   { 
      rinexFile.setPathFilenameMode( rinexFilename, ios::in ); 
   }
   catch( RinexFileException &openExcep )
   {
      cout << "Error opening the file: " << rinexFilename << endl 
           << "Rinex File Exception: " << endl << openExcep.getMessage() << endl << endl;
   }

   // get the file type
   return rinexFile.getRinexFileType();
}


bool NRinexUtils::OpenRinexFileForInput( NGSrinex::RinexFile &rinexFile, const std::string& rinexFilename )
{
   // open the input RINEX file
   try
   {
      rinexFile.setPathFilenameMode( rinexFilename, ios_base::in );
   }
   catch( ... )
   {
      return false;
   }

   // success
   return true;
}


bool NRinexUtils::OpenRinexObservationFileForInput( NGSrinex::RinexObsFile &rinexFile, const std::string& rinexFilename )
{
   // open the file
   if( !OpenRinexFileForInput( rinexFile, rinexFilename ) )
      return false;

   // read the header
   try
   {
      rinexFile.readHeader();
   }
   catch( ... )
   {
      return false;
   }

   // return success
   return true;
}


bool NRinexUtils::OpenRinexNavigationFileForInput( NGSrinex::RinexNavFile &rinexFile, const std::string& rinexFilename )
{
   // open the file
   if( !OpenRinexFileForInput( rinexFile, rinexFilename ) )
      return false;

   // read the header
   try
   {
      rinexFile.readHeader();
   }
   catch( ... )
   {
      return false;
   }

   // return success
   return true;
}


bool NRinexUtils::OpenRinexFileForOutput( NGSrinex::RinexFile &rinexFile, const std::string& rinexFilename )
{
   // open the input RINEX file
   try
   {
      rinexFile.setPathFilenameMode( rinexFilename, ios_base::out );
   }
   catch( ... )
   {
      return false;
   }

   // success
   return true;
}


