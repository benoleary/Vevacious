/*
 * VevaciousBatchRunner.cpp
 *
 *  Created on: Mar 12, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of Vevacious, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.Vevacious.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */


#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <dirent.h>
#include "BOLlib/include/FilePlaceholderManager.hpp"
#include "VevaciousRunner.hpp"
#include "SLHA.hpp"


int main( int argumentCount,
          char** argumentCharArrays )
{
  BOL::ArgumentParser argumentParser( argumentCount,
                                      argumentCharArrays,
                                      "input",
                                      "VevaciousInitialization.xml" );

  BOL::FilePlaceholderManager placeholderManager( "",
                                                  ".placeholder",
                                                  "" );

  Vevacious::VevaciousRunner vevaciousRunner( argumentParser );

  std::string mainPythonFilename( argumentParser.fromTag( "python_main",
                                                          "Vevacious.py" ) );
  std::string resultFilename( argumentParser.fromTag( "result_file",
                                                      "Vevacious.vout" ) );

  std::string inputDirectory(  argumentParser.fromTag( "input_dir",
                                                       "./" ) );
  std::string outputDirectory(  argumentParser.fromTag( "output_dir",
                                                        "./" ) );

  if( outputDirectory.compare( inputDirectory ) == 0 )
  {
    std::cout
    << std::endl
    << "Input directory and output directory must be different, as the"
    << " filenames do not change!";
    std::cout << std::endl;

    return EXIT_FAILURE;
  }

  placeholderManager.prepareFilenames( inputDirectory,
                                       outputDirectory,
                                       outputDirectory );

  std::ifstream inputSlha;
  std::string slhaLine( "" );
  std::ofstream outputSlha;

  while( placeholderManager.holdNextPlace() )
  {
    std::cout
    << std::endl
    << "Looking for input file \"" << placeholderManager.getInput()
    << "\", using placeholder \"" << placeholderManager.getPlaceholder()
    << "\", with result to go to \"" << placeholderManager.getOutput() << "\"";
    std::cout << std::endl;

    try
    {
      vevaciousRunner.findTreeLevelExtrema( placeholderManager.getInput() );
      vevaciousRunner.prepareParameterDependentPython(
                                               placeholderManager.getInput() );
    }
    catch( std::invalid_argument& invalidSlha )
    {
      std::cout
      << std::endl
      << "Unfortunately the SLHA file \"" << placeholderManager.getInput()
      << "\" was not acceptable (" << invalidSlha.what()
      << "), hence the calculation has been aborted.";
      std::cout << std::endl;

      continue;
    }
    vevaciousRunner.runPython( mainPythonFilename );
    if( !(BOL::UsefulStuff::fileExists( resultFilename )) )
    {
      std::cout
      << std::endl
      << "No result from this SLHA file (\"" << placeholderManager.getInput()
      << "\")!";
      std::cout << std::endl;
      continue;
    }
    inputSlha.open( placeholderManager.getInput().c_str() );
    outputSlha.open( placeholderManager.getOutput().c_str() );
    while( inputSlha.good() )
    {
      std::getline( inputSlha,
                    slhaLine );
      if( slhaLine.find( "VEVACIOUSRESULT" ) != std::string::npos )
      {
        break;
      }
      outputSlha << slhaLine << std::endl;
    }
    inputSlha.close();
    outputSlha.close();
    vevaciousRunner.appendResultsToSlha( placeholderManager.getOutput() );

    std::cout
    << std::endl
    << "Successfully wrote result in \"" << placeholderManager.getOutput()
    << "\"";
    std::cout << std::endl;
  }

  std::cout
  << std::endl
  << "Finished.";
  std::cout << std::endl;

  // this was a triumph! I'm making a note here:
  return EXIT_SUCCESS;
}

