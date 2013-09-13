/*
 * TadpoleSolver.cpp
 *
 *  Created on: Oct 2, 2012
 *      Authors: José Eliel Camargo (elielx@gmail.com),
 *               Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 José Eliel Camargo, Ben O'Leary
 *
 *      This file is part of Vevacious, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.Vevacious.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "../include/TadpoleSolver.hpp"

namespace Vevacious
{
  std::string const
  TadpoleSolver::relativeHom4ps2InputFilename( "current_tadpoles.in" );
  std::string const
  TadpoleSolver::relativeHom4ps2OutputFilename( "data.roots" );

  TadpoleSolver::TadpoleSolver( std::string const& hom4ps2Directory,
                                std::string const& homotopyType ) :
    hom4ps2Directory( hom4ps2Directory + "/" ),
    homotopyType( homotopyType ),
    absoluteHom4ps2InputFilename( hom4ps2Directory
                                  + relativeHom4ps2InputFilename ),
    absoluteHom4ps2OutputFilename( hom4ps2Directory
                                   + relativeHom4ps2OutputFilename ),
    parsingStream(),
    tadpoleEquations(),
    currentConvertedTadpole(),
    convertedTadpoles(),
    tadpoleSolutionsFile( "###",
                          false ),
    complexSolutions(),
    variableNames(),
    purelyRealSolutionSets(),
    purelyRealSolutionSetCount( 0 ),
    notPurelyRealSolutionSets(),
    notPurelyRealSolutionSetCount( 0 )
  {
    parsingStream.precision( 17 );
  }

  TadpoleSolver::~TadpoleSolver()
  {
    // does nothing.
  }


  void TadpoleSolver::writeTreeLevelSolutions(
                                            std::string const& outputFilename )
  // this writes the solutions from HOM4PS2 in a Mathematica-friendly format
  // in outputFilename.
  {
    std::ofstream outputFile( BOL::StringParser::ensureDirectoryExists(
                                                    outputFilename ).c_str() );
    if( !(outputFile.good()) )
    {
      std::cout
      << std::endl
      << "error! could not open \"" << outputFilename
      << "\" for writing tree-level extrema!";
      std::cout << std::endl;
      throw std::runtime_error( "could not write to " + outputFilename );
    }
    unsigned int const numberOfVariables( variableNames.size() );
    std::stringstream parsingStream;
    parsingStream.precision( 17 );
    std::string parsingString( "" );
    outputFile
    << "Vevacious_NumberOfRealSolutionSetsFoundByHom4ps2IncludingDuplicates = "
    << purelyRealSolutionSetCount
    << "\nVevacious_TreeLevelRealExtrema = { ";
    for( std::set< std::vector< long double > >::iterator
         whichSolutionSet( purelyRealSolutionSets.begin() );
         purelyRealSolutionSets.end() != whichSolutionSet;
         ++whichSolutionSet )
    {
      if( purelyRealSolutionSets.begin() != whichSolutionSet )
      {
        outputFile << ", ";
      }
      outputFile
      << "{ ";
      for( unsigned whichVariable( 0 );
           numberOfVariables > whichVariable;
           ++whichVariable )
      {
        if( 0 < whichVariable )
        {
          outputFile << ", ";
        }
        outputFile
        << variableNames[ whichVariable ] << " -> "
        << longDoubleForMathematica( (*whichSolutionSet)[ whichVariable ] );
      }
      outputFile
      << " }";
    }
    outputFile
    << " }" << std::endl << std::endl;

    outputFile
    << "Vevacious_NumberOfComplexRootsFoundByHom4ps2IncludingDuplicates = "
    << notPurelyRealSolutionSetCount
    << "\nVevacious_TreeLevelComplexExtrema = { ";
    for( std::set< std::vector< std::pair< long double,
                                           long double > > >::iterator
         whichSolutionSet( notPurelyRealSolutionSets.begin() );
         notPurelyRealSolutionSets.end() != whichSolutionSet;
         ++whichSolutionSet )
    {
      if( notPurelyRealSolutionSets.begin() != whichSolutionSet )
      {
        outputFile << ", ";
      }
      outputFile
      << "{ ";
      for( unsigned whichVariable( 0 );
           numberOfVariables > whichVariable;
           ++whichVariable )
      {
        if( 0 < whichVariable )
        {
          outputFile << ", ";
        }
        outputFile
        << variableNames[ whichVariable ] << " -> ( "
        << longDoubleForMathematica(
                                   (*whichSolutionSet)[ whichVariable ].first )
        << " + ( I * " << longDoubleForMathematica(
                       (*whichSolutionSet)[ whichVariable ].second ) << " ) )";
      }
      outputFile
      << " }";
    }
    outputFile
    << " }";
    outputFile.close();
  }

  void TadpoleSolver::writeSpecificTadpoles(
                                           SarahInterpreter& sarahInterpreter )
  // this creates an input file for HOM4PS2 using tadpoleEquations &
  // sarahInterpreter.
  {
    std::ofstream hom4ps2InputFile( absoluteHom4ps2InputFilename.c_str() );
    if( !(hom4ps2InputFile.good()) )
    {
      std::cout
      << std::endl
      << "Error! was not able to write \"" << absoluteHom4ps2InputFilename
      << "\" as the temporary input file for HOM4PS2!";
      std::cout << std::endl;
      throw std::runtime_error( "could not write to "
                                + absoluteHom4ps2InputFilename );
    }
    convertedTadpoles.clear();
    for( int whichTadpole( 0 );
         tadpoleEquations.getLastIndex() > whichTadpole;
         ++whichTadpole )
      // this loop doesn't look at the last string of tadpoleEquations, because
      // it doesn't end in ';'.
    {
      convertTadpoleEquation( BOL::StringParser::trimFromFrontAndBack(
                                              tadpoleEquations[ whichTadpole ],
                                BOL::StringParser::whitespaceAndNewlineChars ),
                              sarahInterpreter );
      if( !(currentConvertedTadpole.empty()) )
      {
        convertedTadpoles.push_back( currentConvertedTadpole );
      }
    }
    if( convertedTadpoles.empty() )
    {
      std::cout
      << std::endl
      << "Error! no tadpole equations (that are not \"0 = 0\")!";
      std::cout << std::endl;
      throw std::runtime_error( "no valid tadpole equations" );
    }

    hom4ps2InputFile << "{ ";
    for( unsigned int whichTadpole( 0 );
         convertedTadpoles.size() > whichTadpole;
         ++whichTadpole )
      // this loop doesn't look at the last string of tadpoleEquations, because
      // it doesn't end in ';'.
    {
      hom4ps2InputFile
      << convertedTadpoles[ whichTadpole ] << ";";
      if( ( convertedTadpoles.size() - 1 ) == whichTadpole )
      {
        hom4ps2InputFile << " }";
      }
      hom4ps2InputFile << std::endl;
    }
    hom4ps2InputFile.close();
  }

  void TadpoleSolver::runHomotopy()
  // this makes a note of the current working directory in
  // originalWorkingDirectory, then changes directory to hom4ps2Directory,
  // then runs HOM4PS2 with relativeHom4ps2InputFilename as input (because it
  // is already in the HOM4PS2 directory), with homotopyType as well, then
  // returns to originalWorkingDirectory.
  {
    char originalWorkingDirectory[ PATH_MAX ];
    if( NULL == getcwd( originalWorkingDirectory,
                        PATH_MAX ) )
    {
      std::cout
      << std::endl
      << "Error! unable to determine current working directory! (necessary,"
      << " since this program needs to change directory to the directory where"
      << " the hom4ps2 executable is, since unfortunately HOM4PS2 runs with"
      << " relative paths; this program returns to where it was called though,"
      << " to make batch calls easier.)";
      std::cout << std::endl;
      throw std::runtime_error(
                             "could not determine current working directory" );
    }
    int directoryChangeSuccess( chdir( hom4ps2Directory.c_str() ) );
    if( 0 != directoryChangeSuccess )
    {
      throw std::runtime_error(
                           "could not change directory to HOM4PS2 directory" );
    }
    std::string systemCommand( "rm ./bin/input.num" );
    BOL::UsefulStuff::runSystemCommand( systemCommand );
    systemCommand.assign( "/bin/bash -c \"./hom4ps2 "
                          + relativeHom4ps2InputFilename + " <<< "
                          + homotopyType + "\"" );
    BOL::UsefulStuff::runSystemCommand( systemCommand );
    // at this point, we are in the directory with hom4ps2 & data.roots.
    directoryChangeSuccess = chdir( originalWorkingDirectory );
    if( 0 != directoryChangeSuccess )
    {
      throw std::runtime_error(
                      "could not change directory back to initial directory" );
    }
  }

  void TadpoleSolver::parseSolutions( SarahInterpreter& sarahInterpreter,
                                      double const zeroImaginaryTolerance )
  // this parses the output of HOM4PS2 into 2 set of sets of solutions of the
  // system of tadpole equations: 1 set where each set has only real
  // solutions, & the other set has the rest of the sets. the solution sets
  // are stored as long doubles or pairs of long doubles.
  {
    complexSolutions.clear();
    variableNames.clear();
    purelyRealSolutionSets.clear();
    notPurelyRealSolutionSets.clear();
    bool successfulOperation( tadpoleSolutionsFile.openFile(
                                             absoluteHom4ps2OutputFilename ) );
    if( !successfulOperation )
    {
      throw std::runtime_error( "could not write to "
                                + absoluteHom4ps2OutputFilename );
    }
    // 1st we pick out lines corresponding to solutions until we find
    // "The order of variables :" which comes after all solutions have been
    // printed.
    std::string lineString( "" );
    std::pair< long double, long double > currentComplexNumber( 0.0L,
                                                                0.0L );
    while( tadpoleSolutionsFile.readNextNonEmptyLineOfFileWithoutComment(
                                                                 lineString ) )
    {
      if( '(' == lineString[ 0 ] )
      {
        // if the above conditions are satisfied, lineString now contains the
        // root as a complex number, in the form where zero is
        // "(  0.0000000000000000E+00 ,  0.0000000000000000E+00)"
        BOL::StringParser::substituteCharacterWith( lineString,
                                                    ',',
                                                    ' ' );
        parsingStream.clear();
        parsingStream.str( lineString.substr( 1,
                                              ( lineString.size() - 2 ) ) );
        parsingStream
        >> currentComplexNumber.first >> currentComplexNumber.second;
        complexSolutions.push_back( currentComplexNumber );
      }
      else if( 0 == lineString.compare( "The order of variables :" ) )
      {
        break;
      }
    }
    // at this point the line "The order of variables :" should have been
    // found. if it hasn't, the file is malformed, but we carry on regardless,
    // looking for the variables in order:
    while( tadpoleSolutionsFile.readNextNonEmptyLineOfFileWithoutComment(
                                                                 lineString ) )
    {
      if( 0 == lineString.compare(
                         "===============>   HOM4PS-2.0   <===============" ) )
      {
        break;
      }
      if( 1 != lineString.size() )
      {
        std::cout
        << std::endl
        << "WARNING! HOM4PS-2.0 reported a name for a VEV that is not a single"
        << " character, which should not have happened, as all VEV names"
        << " should have been converted to single-character names before"
        << " HOM4PS-2.0 was run! The offending name in question is \""
        << lineString << "\".";
        std::cout << std::endl;
      }
      variableNames.push_back( lineString[ 0 ] );
    }
    // at this point, variableNames has all the names of the variables
    // appearing between the lines "The order of variables :" &
    // "===============>   HOM4PS-2.0   <===============".
    sarahInterpreter.sortVevsToMatchSolutions( variableNames );
    unsigned int const numberOfVariables( variableNames.size() );

    // now we read in the numbers of roots as a check:
    tadpoleSolutionsFile.readNextNonEmptyLineOfFileWithoutComment(
                                                                  lineString );
    // now lineString should be "The # of roots      = " + the total number of
    // roots.
    parsingStream.clear();
    parsingStream.str( lineString.substr( ( lineString.find( '=' )
                                            + 1 ) ) );
    unsigned int totalNumberOfRoots;
    parsingStream >> totalNumberOfRoots;
    if( complexSolutions.size() != ( numberOfVariables * totalNumberOfRoots ) )
    {
      std::cout
      << std::endl
      << "warning! read in " << complexSolutions.size()
      << " complex numbers, should have read in " << variableNames.size()
      << " (the number of variables) * " << totalNumberOfRoots
      << " (the total number of solutions declared by HOM4PS2) = "
      << ( variableNames.size() * totalNumberOfRoots ) << " complex numbers!";
      std::cout << std::endl;
    }

    tadpoleSolutionsFile.readJustNextValidLine( lineString );
    // now lineString should be "The # of real roots      = " + the number of
    // real roots.
    parsingStream.clear();
    parsingStream.str( lineString.substr( ( lineString.find( '=' )
                                            + 1 ) ) );
    int numberOfRealRoots;
    parsingStream >> numberOfRealRoots;
    tadpoleSolutionsFile.closeFile();

    currentComplexNumber.first = BOL::UsefulStuff::notANumber;
    currentComplexNumber.second = BOL::UsefulStuff::notANumber;
    std::vector< long double > currentRealSolution( numberOfVariables,
                                                BOL::UsefulStuff::notANumber );
    std::vector< std::pair< long double, long double > >
    currentComplexSolution( numberOfVariables,
                            currentComplexNumber );
    unsigned int whichVariable( 0 );
    bool setWasPurelyReal( true );
    purelyRealSolutionSetCount = 0;
    notPurelyRealSolutionSetCount = 0;
    for( std::vector< std::pair< long double, long double > >::iterator
         whichRoot( complexSolutions.begin() );
         complexSolutions.end() > whichRoot;
         ++whichRoot )
    {
      currentRealSolution[ whichVariable ] = whichRoot->first;
      currentComplexSolution[ whichVariable ] = *whichRoot;
      if( zeroImaginaryTolerance < abs( whichRoot->second ) )
      {
        setWasPurelyReal = false;
      }
      ++whichVariable;

      if( !( numberOfVariables > whichVariable ) )
      {
        // when we've read in a full set of solutions, we note either a purely
        // real solution set or a solution set that was not purely real:
        if( setWasPurelyReal )
        {
          purelyRealSolutionSets.insert( currentRealSolution );
          ++purelyRealSolutionSetCount;
        }
        else
        {
          notPurelyRealSolutionSets.insert( currentComplexSolution );
          ++notPurelyRealSolutionSetCount;
        }
        // either way, we reset setWasPurelyReal for the next solution set &
        // note that we are starting a new set of solutions:
        whichVariable = 0;
        setWasPurelyReal = true;
      }
    }

    if( purelyRealSolutionSetCount != numberOfRealRoots )
    {
      std::cout
      << std::endl
      << "warning! read in " << purelyRealSolutionSetCount
      << " sets of complex numbers all with zero imaginary part, HOM4PS2"
      << " declared " << numberOfRealRoots << " sets of real numbers!";
      std::cout << std::endl;
    }
  }

  void TadpoleSolver::convertTadpoleEquation( std::string const& tadpoleString,
                                           SarahInterpreter& sarahInterpreter )
  // this puts tadpoleString with its SLHA values substituted in into
  // currentConvertedTadpole.
  {
    size_t termStart( tadpoleString.find( "(" ) );
    size_t termEnd( tadpoleString.find( ")",
                                        termStart ) );
    currentConvertedTadpole.assign( "" );
    std::string convertedTerm( "" );
    while( std::string::npos != termEnd )
    {
      convertedTerm.assign( sarahInterpreter( tadpoleString.substr(
                                                             ( termStart + 1 ),
                                                 ( termEnd - termStart - 1 ) ),
                                              true,
                                              true ) );
      if( !(convertedTerm.empty()) )
      {
        if( !(currentConvertedTadpole.empty()) )
        {
          currentConvertedTadpole.append( " + " );
        }
        currentConvertedTadpole.append( "(" + convertedTerm + ")" );
      }
      termStart = tadpoleString.find( "(",
                                      termEnd );
      termEnd = tadpoleString.find( ")",
                                    termStart );
    }
  }

} /* namespace Vevacious */
