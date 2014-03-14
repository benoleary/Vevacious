/*
 * Vevacious.cpp
 *
 *  Created on: Oct 2, 2012
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
#include <ctime>
#include <stdexcept>
#include <sys/time.h>
#include "BOLlib/include/ArgumentParser.hpp"
#include "BOLlib/include/StringParser.hpp"
#include "BOLlib/include/UsefulStuff.hpp"
#include "VevaciousRunner.hpp"

double secondsSince( timeval const& referenceTimeval )
{
  timeval currentTimeval;
  gettimeofday( &currentTimeval,
                NULL );
  return ( (double)( currentTimeval.tv_sec - referenceTimeval.tv_sec )
           + 0.000001
             * (double)( currentTimeval.tv_usec - referenceTimeval.tv_usec ) );
}


int main( int argumentCount,
          char** argumentCharArrays )
{
  // the timing of each part of the calculation is displayed in case anyone is
  // curious.
  timeval startTimeval;
  gettimeofday( &startTimeval,
                NULL );

  BOL::ArgumentParser argumentParser( argumentCount,
                                      argumentCharArrays,
                                      "input",
                                      "VevaciousInitialization.xml" );

  // set up runner for specific model
  Vevacious::VevaciousRunner
  vevaciousRunner( argumentParser );

  // this setting almost certainly has to be changed from the default ( "./" ):
  vevaciousRunner.setPathToCosmotransitions( argumentParser.fromTag( "ct_path",
                                                     "./CosmoTransitions/" ) );

  // the other settings have "reasonable" defaults, but one might want to
  // change them here:
  std::string resultsFilename( argumentParser.fromTag( "result_file",
                                                       "./MyResult.vout" ) );
  vevaciousRunner.setResultsFilename( resultsFilename );
  BOL::UsefulStuff::runSystemCommand( "rm " + resultsFilename );
  std::string
  maxSaddleNudgesInput( argumentParser.fromTag( "max_saddle_nudges",
                                                "" ) );
  if( !(maxSaddleNudgesInput.empty()) )
  {
    vevaciousRunner.setMaximumMinuitNudgesOffSaddlePoints(
                                                        maxSaddleNudgesInput );
  }

  double setupFinishTime( secondsSince( startTimeval ) );
  std::cout
  << std::endl
  <<
  "   @@@...@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n"
  "   @@@@@@.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n"
  "   @@@@@@?.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n"
  "   @@@@@@@.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n"
  "   @@@@@@@..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n"
  "   @@@@@@@..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n"
  "   @@@@@@@~ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n"
  "   @@@@@@@?.?@@@@@@@@@@@@@@@@@@@@@@.~@@@@ \n"
  "   @@@@@@@@..@@@@@@@@@@@@@@@@@@@@..@@@@@@ \n"
  "   @@@@@@@@..@@@@@@@@@@@@@@@@@@@,.@@@@@@@ \n"
  "   @@@@@@@@..@@@@@@@@@@@@@@@@@@@.+@@@@@@@ \n"
  "   @@@@@@@@..@@@@@@@@@@@@@@@@@@@ @@@@@@@@ \n"
  "   @@@@@@@@ .@@@@@@@@@@@@@@@@@@. @@@@@@@@ \n"
  "   @@@@@@@@..~@@@@@@@@@@@@@@@@@..@@@@@@@@ \n"
  "   @@@@@@@@.  @@@@@@@@@@@@@@@@@..@@@@@@@@ \n"
  "   @@@@@@@@. .@@@@@@@@@@@@@@@@@..@@@@@@@@ \n"
  "   @@@@@@@@.  ~@@@@@@@@@@@@@@@. .@@@@@@@@ \n"
  "   @@@@@@@@?. .@@@@@@@@@@@@@@@. @@@@@@@@@ \n"
  "   @@@@@@@@@.  .@@@@@@@@@@@@@.  @@@@@@@@@ \n"
  "   @@@@@@@@@.. .@@@@@@@@@@@@@. .@@@@@@@@@ \n"
  "   @@@@@@@@@@.  .@@@@###@@@@.  @@@@@@@@@@ \n"
  "   @@@@@@@@@@@.  .@@#####@@....@@@@@@@@@@ \n"
  "   @@@@@@@@@@@,.  ..@###@@....@@@@@@@@@@@ \n"
  "   @@@@@@@@@@@@,.  ..,@@.   .@@@@@@@@@@@@ \n"
  "__   _______   ____ _  ___ _  ___  _   _ ___ \n"
  "\\ \\ / / _ \\ \\ / / _` |/ __| |/ _ \\| | | / __|\n"
  " \\ V /  __/\\ V / (_| | (__| | (_) | |_| \\__ \\\n"
  "  \\_/ \\___| \\_/ \\__,_|\\___|_|\\___/ \\__,_|___/\n"
  "                                             \n"
  "\n"
  "\n"
  "\n"
  << "Setting up the VevaciousRunner object took " << setupFinishTime
  << " seconds.";
  std::cout << std::endl;

  // the initial points for PyMinuit do not _have_ to be the extrema of the
  // tree-level potential of the Lagrangian density used by PyMinuit! it's
  // probably best to consistently use only 1 SLHA file for both parts, but
  // one might have a motivation for using 2 SLHA files: say file A, which has
  // the values of the parameters that one wants to base the one-loop results
  // on, and file B, which leads to a tree-level potential with extrema close
  // to where one thinks the one-loop potential for file A has its minima.
  // anyway, Vevacious allows one to use different SLHA files for each part.

  // solve tadpoles for specific SLHA, adding results to given (/default) file.
  std::string slhaFilename( argumentParser.fromTag( "slha_file",
                                                 "./MyParameters.slha.out" ) );
  try
  {
    vevaciousRunner.findTreeLevelExtrema( slhaFilename );
  }
  catch( std::invalid_argument& invalidSlha )
  {
    std::cout
    << std::endl
    << "Unfortunately the SLHA file was not acceptable ("
    << invalidSlha.what() << "), hence the calculation has been aborted.";
    std::cout << std::endl;

    return EXIT_FAILURE;
  }

  // this step is not necessary, but the tree-level extrema are printed out for
  // reference in a Mathematica-friendly format.
  BOL::UsefulStuff::runSystemCommand(
                                     "rm ./Vevacious_tree-level_extrema.txt" );

  vevaciousRunner.writeCalculatedTreeLevelExtrema(
                                        "./Vevacious_tree-level_extrema.txt" );

  double extremaFinishTime( secondsSince( startTimeval ) );
  std::cout
  << std::endl
  << "Finding all the tree-level extrema and parsing the results took "
  << ( extremaFinishTime - setupFinishTime ) << " seconds.";
  std::cout << std::endl;

  // one could append more extrema using another SLHA file as well, if there
  // was a reason to think that those extrema might help PyMinuit...

  // prepare potential for Vevacious.py (/given other Python program).
  try
  {
    vevaciousRunner.prepareParameterDependentPython( slhaFilename );
  }
  catch( std::invalid_argument& invalidSlha )
  {
    std::cout
    << std::endl
    << "Unfortunately the SLHA file was not acceptable ("
    << invalidSlha.what() << "), hence the calculation has been aborted.";
    std::cout << std::endl;

    return EXIT_FAILURE;
  }

  // run Vevacious.py (/given other Python program).
  std::string mainPythonFilename( argumentParser.fromTag( "python_main",
                                                          "./Vevacious.py" ) );

  vevaciousRunner.runPython( mainPythonFilename );

  double minimizingFinishTime( secondsSince( startTimeval ) );
  std::cout
  << std::endl
  << "Minimizing the one-loop potential, calculating the tunneling time or its"
  << " upper bound, and printing the results took "
  << ( minimizingFinishTime - extremaFinishTime ) << " seconds.";
  std::cout << std::endl;

  // append the results to the SLHA file in custom blocks:
  std::string howToAppendSlhaBlock( argumentParser.fromTag( "appendToSlha",
                                                            "normal" ) );
  if( !(howToAppendSlhaBlock.empty()) )
  {
    vevaciousRunner.appendResultsToSlha( slhaFilename,
                              ( 0 != howToAppendSlhaBlock.compare( "SSP" ) ) );
  }

  std::cout
  << std::endl
  << "Vevacious finished. Total running time was "
  << secondsSince( startTimeval ) << " seconds.";
  std::cout << std::endl;

  // this was a triumph! I'm making a note here:
  return EXIT_SUCCESS;
}


