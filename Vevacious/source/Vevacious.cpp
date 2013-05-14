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
  vevaciousRunner( argumentParser.fromTag( "model_file",
                                          "./Vevacious.in.realStauVevs_MSSM" ),
                   argumentParser.fromTag( "hom4ps2_dir",
                                           "./HOM4PS2/" ),
                   argumentParser.fromTag( "homotopy_type",
                                           "1" ) );

  // this setting almost certainly has to be changed from the default ( "./" ):
  vevaciousRunner.setPathToCosmotransitions( argumentParser.fromTag( "ct_path",
                                                     "./CosmoTransitions/" ) );

  // the other settings have "reasonable" defaults, but one might want to
  // change them here:
  std::string resultsFilename( argumentParser.fromTag( "result_file",
                                "./VevaciousResults.xml.realStauVevs_MSSM" ) );
  vevaciousRunner.setResultsFilename( resultsFilename );
  BOL::UsefulStuff::runSystemCommand( "rm " + resultsFilename );
  vevaciousRunner.setImaginaryPartTolerance( argumentParser.fromTag(
                                                         "imaginary_tolerance",
                                                               "0.0000001" ) );
  vevaciousRunner.setMinuitNudgesOffSaddlePoints(
                                       argumentParser.fromTag( "saddle_nudges",
                                                               "1.0, 2.0" ) );
  vevaciousRunner.setMaximumMinuitNudgesOffSaddlePoints(
                                   argumentParser.fromTag( "max_saddle_nudges",
                                                         "2" ) );
  vevaciousRunner.setMinuitRollingTolerance( argumentParser.fromTag(
                                                              "roll_tolerance",
                                                                     "0.1" ) );
  vevaciousRunner.setLifetimeForDirectPath( argumentParser.fromTag(
                                                                 "direct_time",
                                                                     "0.1" ) );
  vevaciousRunner.setLifetimeForDeformedPath( argumentParser.fromTag(
                                                               "deformed_time",
                                                                     "0.1" ) );
  // negative arguments (as doubles or std::strings representing doubles) for
  // setLifetimeForDirectPath, or setLifetimeForDeformedPath means that the
  // relevant tunneling time calculation will be skipped.

  double setupFinishTime( secondsSince( startTimeval ) );
  std::cout
  << std::endl
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
                                                    "./SPheno.spc.MSSM" ) );
  vevaciousRunner.overwriteTreeLevelExtrema( slhaFilename );

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
  vevaciousRunner.prepareParameterDependentPython( slhaFilename );

  // run Vevacious.py (/given other Python program).
  std::string mainPythonFilename( argumentParser.fromTag( "python_main",
                                                          "./Vevacious.py" ) );

  // actually, unfortunately this doesn't work because fork()ing seems to lose
  // the PYTHONPATH.
  /**
  int deformationPatience( BOL::StringParser::stringToInt(
                                argumentParser.fromTag( "deformation_patience",
                                                        "-1" ) ) );
  bool finishedOnTime( vevaciousRunner.runPython( mainPythonFilename,
                                                  deformationPatience ) );
  if( !finishedOnTime )
  {
    std::cout
    << std::endl
    << "Warning! Calculation taking too long: killing Python program and"
    << " running again, but without allowing CosmoTransitions to deform the"
    << " path in VEV-space from the false vacuum to the true vacuum.";
    std::cout << std::endl;

    vevaciousRunner.setLifetimeForDeformedPath( -1.0 );
    vevaciousRunner.prepareParameterDependentPython( slhaFilename );
    vevaciousRunner.runPython( mainPythonFilename );
    vevaciousRunner.setLifetimeForDeformedPath( argumentParser.fromTag(
                                                               "deformed_time",
                                                                     "0.1" ) );
  }
  **/

  // ... so we just have to use the version that doesn't try to kill the Python
  // program no matter how long it takes.
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


