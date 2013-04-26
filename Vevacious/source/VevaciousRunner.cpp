/*
 * VevaciousRunner.cpp
 *
 *  Created on: Mar 5, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of Vevacious, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.Vevacious.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "../include/VevaciousRunner.hpp"

namespace Vevacious
{
  double const VevaciousRunner::lifetimeFactor( ( 4.0 * 44.0 * log( 10.0 ) )
                                                + log( 0.0001 ) );
  // The age of the Universe is pretty close to 10^(44)/TeV. The equation used
  // for the tunneling time t is t^(-4) = A exp( -B ), where B is the action &
  // A is a dimensionful factor depending on the solitonic solutions. It is
  // expected to be of the order of (0.1 TeV)^4. For a given t, B = ln( A t^4 )
  // and if t is r * 10^(44)/TeV, then
  // B = 4 ln( r ) + 4 ln( 10^(44) ) + ln( A * TeV^(-4) ), explaining the above
  // form of lifetimeFactor.
  std::string const
  VevaciousRunner::pointsToTry( "pointsToTry" );
  std::string const
  VevaciousRunner::resultsFilenameVariableName( "outputFile" );
  std::string const
  VevaciousRunner::maximumSaddleSplitsVariableName( "maximumSaddleSplits" );
  std::string const
  VevaciousRunner::rollingToleranceVariableName( "rollingTolerance" );
  std::string const
  VevaciousRunner::quarticActionBoundVariableName( "quarticActionBound" );
  std::string const VevaciousRunner::pathToCosmotransitionsVariableName(
                                                    "pathToCosmotransitions" );
  std::string const
  VevaciousRunner::directActionBoundVariableName( "directActionBound" );
  std::string const
  VevaciousRunner::deformedActionBoundVariableName( "deformedActionBound" );
  std::string const VevaciousRunner::defaultPythonFilename( "Vevacious.py" );
  std::string const VevaciousRunner::treeLevelExtremaHeader(
         "import VevaciousParameterDependent as VPD\n\npointsToTry = []\n\n" );
  BOL::StringParser const VevaciousRunner::slhaIndexMaker( 3,
                                                           ' ',
                                                           9,
                                                           3,
                                                           "" );
  BOL::StringParser const VevaciousRunner::slhaDoubleMaker( 9,
                                                            ' ',
                                                            9,
                                                            3 );

  VevaciousRunner::VevaciousRunner( std::string const& modelFilename,
                                    std::string const& hom4ps2Directory,
                                    std::string const& homotopyType ) :
    sarahInterpreter(),
    tadpoleSolver( hom4ps2Directory,
                   homotopyType ),
    potentialMinimizer(),
    treeLevelExtrema(),
    firstWriteOfExtrema( true ),
    resultsFilename( "VevaciousResult.xml" ),
    imaginaryTolerance( 0.000001 ),
    maximumSaddleSplits( 2 ),
    rollingTolerance( 0.2 ),
    quarticActionBound( actionFromLifetime( 0.1 ) ),
    pathToCosmotransitions( "./" ),
    directActionBound( quarticActionBound ),
    deformedActionBound( quarticActionBound )
  {
    BOL::AsciiXmlParser xmlParser;
    bool xmlOpened( xmlParser.openRootElementOfFile( modelFilename ) );
    if( !xmlOpened )
    {
      std::cout
      << std::endl
      << "Error! could not open \"" << modelFilename
      << "\", so there is no potential to minimize.";
      std::cout << std::endl;
    }
    else
    {
      while( xmlParser.readNextElement() )
      {
        if( xmlParser.currentElementNameMatches( "input_vevs" ) )
        {
          sarahInterpreter.setVevNames(
                                     xmlParser.getCurrentElementAttributes() );
          sarahInterpreter.setPositiveVevs(
                                 xmlParser.getTrimmedCurrentElementContent() );
        }
        else if( xmlParser.currentElementNameMatches( "block_prefixes" ) )
        {
          sarahInterpreter.setBlockPrefixes(
                                     xmlParser.getCurrentElementAttributes() );
        }
        else if( xmlParser.currentElementNameMatches( "tadpoles" ) )
        {
          tadpoleSolver.setEquations(
                                 xmlParser.getTrimmedCurrentElementContent() );
        }
        else if( xmlParser.currentElementNameMatches( "polynomial_part" ) )
        {
          potentialMinimizer.setPolynomialPartOfPotential(
                                 xmlParser.getTrimmedCurrentElementContent() );
        }
        else if( xmlParser.currentElementNameMatches( "mass-squared_matrix" ) )
        {
          potentialMinimizer.addMassSquaredMatrix(
                                       xmlParser.getCurrentElementAttributes(),
                                 xmlParser.getTrimmedCurrentElementContent() );
        }
      }
    }
  }

  VevaciousRunner::VevaciousRunner( BOL::ArgumentParser& argumentParser ) :
    sarahInterpreter(),
    tadpoleSolver( argumentParser.fromTag( "hom4ps2_dir",
                                           "./HOM4PS2/" ),
                   argumentParser.fromTag( "homotopy_type",
                                           "1" ) ),
    potentialMinimizer(),
    treeLevelExtrema(),
    firstWriteOfExtrema( true ),
    resultsFilename( argumentParser.fromTag( "result_file",
                                "./VevaciousResults.xml.realStauVevs_MSSM" ) ),
    imaginaryTolerance( BOL::StringParser::stringToDouble(
                                 argumentParser.fromTag( "imaginary_tolerance",
                                                         "0.0000001" ) ) ),
    maximumSaddleSplits( BOL::StringParser::stringToInt(
                                   argumentParser.fromTag( "max_saddle_nudges",
                                                           "3" ) ) ),
    rollingTolerance( BOL::StringParser::stringToDouble(
                                      argumentParser.fromTag( "roll_tolerance",
                                                              "0.1" ) ) ),
    quarticActionBound( actionFromLifetime( BOL::StringParser::stringToDouble(
                                        argumentParser.fromTag( "quartic_time",
                                                                "-0.1" ) ) ) ),
    pathToCosmotransitions( argumentParser.fromTag( "ct_path",
                                                    "./CosmoTransitions/" ) ),
    directActionBound( actionFromLifetime( BOL::StringParser::stringToDouble(
                                         argumentParser.fromTag( "direct_time",
                                                                 "0.1" ) ) ) ),
    deformedActionBound( actionFromLifetime( BOL::StringParser::stringToDouble(
                                       argumentParser.fromTag( "deformed_time",
                                                               "0.1" ) ) ) )
  {
    BOL::AsciiXmlParser xmlParser;
    std::string modelFilename( argumentParser.fromTag( "model_file",
                                        "./Vevacious.in.realStauVevs_MSSM" ) );
    bool xmlOpened( xmlParser.openRootElementOfFile( modelFilename ) );
    if( !xmlOpened )
    {
      std::cout
      << std::endl
      << "Error! could not open \"" << modelFilename
      << "\", so there is no potential to minimize.";
      std::cout << std::endl;
    }
    else
    {
      while( xmlParser.readNextElement() )
      {
        if( xmlParser.currentElementNameMatches( "input_vevs" ) )
        {
          sarahInterpreter.setVevNames(
                                     xmlParser.getCurrentElementAttributes() );
          sarahInterpreter.setPositiveVevs(
                                 xmlParser.getTrimmedCurrentElementContent() );
        }
        else if( xmlParser.currentElementNameMatches( "block_prefixes" ) )
        {
          sarahInterpreter.setBlockPrefixes(
                                     xmlParser.getCurrentElementAttributes() );
        }
        else if( xmlParser.currentElementNameMatches( "tadpoles" ) )
        {
          tadpoleSolver.setEquations(
                                 xmlParser.getTrimmedCurrentElementContent() );
        }
        else if( xmlParser.currentElementNameMatches( "polynomial_part" ) )
        {
          potentialMinimizer.setPolynomialPartOfPotential(
                                 xmlParser.getTrimmedCurrentElementContent() );
        }
        else if( xmlParser.currentElementNameMatches( "mass-squared_matrix" ) )
        {
          potentialMinimizer.addMassSquaredMatrix(
                                       xmlParser.getCurrentElementAttributes(),
                                 xmlParser.getTrimmedCurrentElementContent() );
        }
      }
    }
  }

  VevaciousRunner::~VevaciousRunner()
  {
    // does nothing.
  }


  void VevaciousRunner::writeCalculatedTreeLevelExtrema(
                                            std::string const& outputFilename )
  {
    std::set< std::vector< long double > > const&
    solutionSets( tadpoleSolver.getPurelyRealSolutionSets() );
    std::vector< char > const&
    vevsOrderedBySolution( sarahInterpreter.getVevsOrderedBySolutions() );
    std::ofstream outputFile( outputFilename.c_str() );
    if( !(outputFile.good()) )
    {
      throw std::runtime_error( "could not open \"" + outputFilename
                                + "\" to write tree-level extrema." );
    }
    outputFile << "VevaciousTreeLevelExtrema = {\n";
    for( std::set< std::vector< long double > >::const_iterator
         whichSet( solutionSets.begin() );
         solutionSets.end() != whichSet;
         ++whichSet )
    {
      if( solutionSets.begin() != whichSet )
      {
        outputFile << ", ";
      }
      outputFile << "{ ";
      for( unsigned int whichField( 0 );
           whichSet->size() > whichField;
           ++whichField )
      {
        if( 0 != whichField )
        {
          outputFile << ", ";
        }
        outputFile
        << sarahInterpreter.getUserVevName(
                                          vevsOrderedBySolution[ whichField ] )
        << " -> ( " << (*whichSet)[ whichField ] << " )";
      }
      outputFile << " }\n";
    }
    outputFile << "}\n";
    outputFile.close();
  }

  void VevaciousRunner::prepareParameterDependentPython(
                                               std::string const& slhaFilename,
                                             std::string const pythonFilename )
  {
    sarahInterpreter.setSlhaFile( slhaFilename );
    BOL::UsefulStuff::runSystemCommand( "rm " + pythonFilename + "c" );
    std::ofstream outputFile( pythonFilename.c_str() );
    if( !(outputFile.good()) )
    {
      throw std::runtime_error( "could not open \"" + pythonFilename
                                + "\" to write potential energy functions." );
    }
    outputFile <<
"from __future__ import division\n"
"import math\n"
"import numpy.linalg\n"
<< resultsFilenameVariableName << " = \"" << resultsFilename << "\"\n"
<< maximumSaddleSplitsVariableName << " = " << maximumSaddleSplits << "\n"
<< rollingToleranceVariableName << " = " << rollingTolerance << "\n"
<< quarticActionBoundVariableName << " = " << quarticActionBound << "\n"
<< pathToCosmotransitionsVariableName << " = \"" << pathToCosmotransitions
<< "\"\n"
<< directActionBoundVariableName << " = " << directActionBound << "\n"
<< deformedActionBoundVariableName << " = " << deformedActionBound << "\n"
<< potentialMinimizer.prepareParameterDependentPython( sarahInterpreter );
    outputFile.close();
  }

  void VevaciousRunner::appendResultsToSlha( std::string const& slhaFilename,
                          bool const properFormatRatherThanSspReadable ) const
  // this takes the XML-format results from the file called resultsFilename,
  // & appends the results & warnings in custom SLHA blocks to the end of the
  // file name slhaFilename.
  {
    double stabiltyResult( -2.0 );
    // -2 will do as an error code for the SLHA block for the moment.
    std::string stabilityVerdict( "error" );
    std::map< std::string, std::string > inputMinimumDepthAndVevs;
    std::map< std::string, std::string > globalMinimumDepthAndVevs;
    std::string lifetimeBound( "-2.0" );
    std::string actionCalculation( "error" );
    std::string warningLine( "" );
    std::string warningBlock( "" );
    int warningNumber( 0 );
    BOL::AsciiXmlParser resultParser( false );
    if( resultParser.openRootElementOfFile( resultsFilename ) )
    {
      while( resultParser.readNextElement() )
      {
        if( resultParser.currentElementNameMatches( "stability" ) )
        {
          stabilityVerdict.assign(
                              resultParser.getTrimmedCurrentElementContent() );
          if( 0 == stabilityVerdict.compare( "stable" ) )
          {
            stabiltyResult = 1.0;
          }
          else if( 0 == stabilityVerdict.compare( "metastable" ) )
          {
            stabiltyResult = 0.0;
          }
          else if( 0 == stabilityVerdict.compare( "unstable" ) )
          {
            stabiltyResult = -1.0;
          }
        }
        else if( resultParser.currentElementNameMatches( "global_minimum" ) )
        {
          globalMinimumDepthAndVevs
          = resultParser.getCurrentElementAttributes();
        }
        else if( resultParser.currentElementNameMatches( "input_minimum" ) )
        {
          inputMinimumDepthAndVevs
          = resultParser.getCurrentElementAttributes();
        }
        else if( resultParser.currentElementNameMatches( "lifetime" ) )
        {
          lifetimeBound = slhaDoubleFromQuotedString(
                              resultParser.getTrimmedCurrentElementContent() );
          std::map< std::string, std::string >::const_iterator
          actionFinder( resultParser.getCurrentElementAttributes().find(
                                                      "action_calculation" ) );
          if( resultParser.getCurrentElementAttributes().end()
              != actionFinder )
          {
            actionCalculation.assign( BOL::StringParser::trimFromFrontAndBack(
                                                          actionFinder->second,
                                                                    "\'\"" ) );
          }
        }
        else if( resultParser.currentElementNameMatches( "warning" ) )
        {
          warningLine.assign( resultParser.getTrimmedCurrentElementContent() );
          Vevacious::SarahInterpreter::removeNewlinesFrom( warningLine );
          warningBlock.append( " "
                               + slhaIndexMaker.intToString( ++warningNumber )
                               + "    " + warningLine + "\n" );
        }
      }
      resultParser.closeFile();
      std::fstream outputFile( slhaFilename.c_str() );
      if( !(outputFile.good()) )
      {
        throw std::runtime_error( "could not open \"" + slhaFilename
                                  + "\" to append results." );
      }
      long endPosition( outputFile.seekg( 0,
                                          std::ios::end ).tellg() );
      while( '\n' == (char)(outputFile.seekg( (--endPosition) ).peek()) )
      {
        // this loop just brings the get pointer back to the char before the
        // last '\n'.
      }
      outputFile.seekp( (++endPosition) );
      // the put pointer is now about to overwrite the 1st '\n' of the sequence
      // of '\n' characters ending the file.
      outputFile << "\n"
      << "BLOCK VEVACIOUSRESULTS # results from Vevacious\n"
      << "    0   0    " << slhaDoubleMaker.doubleToString( stabiltyResult );
      if( !properFormatRatherThanSspReadable )
      {
        outputFile << " # ";
      }
      outputFile
      << "    " << stabilityVerdict << "    # stability of input\n"
      << "    0   1    " << lifetimeBound;
      if( !properFormatRatherThanSspReadable )
      {
        outputFile << " # ";
      }
      outputFile << "    " << actionCalculation
      << "    # tunneling time in Universe ages / calculation type\n";
      int vevNumber( -1 );
      for( std::map< std::string, std::string >::const_iterator
           whichVev( inputMinimumDepthAndVevs.begin() );
           inputMinimumDepthAndVevs.end() != whichVev;
           ++whichVev )
      {
        outputFile
        << "    1 " << slhaIndexMaker.intToString( ++vevNumber ) << "    "
        << slhaDoubleFromQuotedString( whichVev->second );
        if( !properFormatRatherThanSspReadable )
        {
          outputFile << " # ";
        }
        outputFile << "    "
        << whichVev->first;
        if( inputMinimumDepthAndVevs.begin() == whichVev )
        {
          outputFile
          << "    # input potential energy density difference from all VEVs ="
          <<                                                          " 0.0\n";
        }
        else
        {
          outputFile
          << "    # input VEV\n";
        }
      }
      vevNumber = -1;
      for( std::map< std::string, std::string >::const_iterator
           whichVev( globalMinimumDepthAndVevs.begin() );
           globalMinimumDepthAndVevs.end() != whichVev;
           ++whichVev )
      {
        outputFile
        << "    2 " << slhaIndexMaker.intToString( ++vevNumber ) << "    "
        << slhaDoubleFromQuotedString( whichVev->second );
        if( !properFormatRatherThanSspReadable )
        {
          outputFile << " # ";
        }
        outputFile << "    "
        << whichVev->first;
        if( globalMinimumDepthAndVevs.begin() == whichVev )
        {
          outputFile
          << "    # global minimum potential energy density difference from"
          <<                                               " all VEVs = 0.0\n";
        }
        else
        {
          outputFile
          << "    # global minimum VEV\n";
        }
      }
      if( properFormatRatherThanSspReadable
          &&
          !(warningBlock.empty()) )
      {
        outputFile
        << "BLOCK VEVACIOUSWARNINGS # warnings from Vevacious\n"
        << warningBlock;
      }
      outputFile.close();
    }
    else
    {
      std::cout
      << std::endl
      << "Warning! \"" << resultsFilename << "\" was not produced!";
      std::cout << std::endl;
    }
  }

  void
  VevaciousRunner::calculateTreeLevelExtrema( std::string const& slhaFilename )
  {
    sarahInterpreter.setSlhaFile( slhaFilename );
    tadpoleSolver.solveTadpoles( sarahInterpreter,
                                 imaginaryTolerance );
    std::set< std::vector< long double > > const&
    solutionSets( tadpoleSolver.getPurelyRealSolutionSets() );
    std::vector< char > const&
    vevsOrderedBySolution( sarahInterpreter.getVevsOrderedBySolutions() );
    treeLevelExtrema.clear();
    treeLevelExtrema.str( "" );
    for( std::set< std::vector< long double > >::const_iterator
         whichSet( solutionSets.begin() );
         solutionSets.end() != whichSet;
         ++whichSet )
    {
      treeLevelExtrema
      << "pointsToTry.append( { ";
      for( unsigned int whichField( 0 );
           whichSet->size() > whichField;
           ++whichField )
      {
        if( 0 != whichField )
        {
          treeLevelExtrema << ", ";
        }
        treeLevelExtrema
        << "\'" << vevsOrderedBySolution[ whichField ] << "\': ( "
        << (*whichSet)[ whichField ] << " / VPD."
        << PotentialMinimizer::energyScale << " )";
      }
      treeLevelExtrema << " } )\n";
    }
    treeLevelExtrema << "\n";
  }

  void VevaciousRunner::writeDefaultPythonProgram() const
  {
    BOL::UsefulStuff::runSystemCommand( "rm " + defaultPythonFilename + "c" );
    std::ofstream outputFile( defaultPythonFilename.c_str() );
    if( !(outputFile.good()) )
    {
      throw std::runtime_error( "could not open \"" + defaultPythonFilename
                                + "\" to write Vevacious Python code." );
    }
    outputFile <<
"# This Python program was automatically generated by Vevacious, a program\n"
"# released under the GNU General Public License, and, as such, this program\n"
"# is also under the GNU General Public License. Please see the\n"
"# README.Vevacious.txt file accompanying the Vevacious C++ code for a full\n"
"# list of files, brief documentation, and further details on the license.\n"
"#\n"
"# Vevacious authors: Ben O'Leary (benjamin.oleary@gmail.com)\n"
"#                    Florian Staub (florian.staub@googlemail.com)\n"
"#                    Jose' Eliel Camargo Molina (elielx@gmail.com)\n"
"#                    Werner Porod (porod@physik.uni-wuerzburg.de)\n"
"#\n"
"#      Copyright 2012 Ben O'Leary, Florian Staub,\n"
"#                     Jose' Eliel Camargo Molina, Werner Porod\n"
"#\n"
"from __future__ import division\n"
"import math\n"
"import numpy\n"
"import numpy.linalg\n"
"import minuit\n"
"import VevaciousParameterDependent as VPD\n"
"import VevaciousTreeLevelExtrema as VTE\n"
"\n"
"warningMessages = []\n"
"namesOfVevs = VPD." << PotentialMinimizer::namesOfVevs << "\n"
"numberOfVevs = len( namesOfVevs )\n"
"stepSize = ( 1.0 / VPD." << PotentialMinimizer::energyScale << " )\n"
"rollingToleranceSquared = ( VPD." << rollingToleranceVariableName << " * VPD."
<< rollingToleranceVariableName << " )\n"
"potentialAtVevOrigin = ( VPD." << PotentialMinimizer::energyScaleFourth
<<                                                                         "\n"
"                         * VPD." << PotentialMinimizer::functionFromDictionary
<<              "( VPD." << PotentialMinimizer::loopCorrectedPotential << ",\n"
"                            VPD." << PotentialMinimizer::vevOrigin << " ) )\n"
"\n"
"# MINUIT\'s hesse() function assumes that it is already at a\n"
"# minimum, but we need to check whether it actually stopped at a\n"
"# saddle point, so we need to work out the Hessian matrix ourselves.\n"
"def NumericalHessian( inputFunction, vevValues ):\n"
"    returnHessian = [ [ 0.0 for vevIndexOne in range( numberOfVevs ) ]\n"
"                                 for vevIndexTwo in range( numberOfVevs ) ]\n"
"    for vevIndexOne in range( numberOfVevs ):\n"
"        firstPoint = vevValues.copy()\n"
"        firstStepSize = stepSize\n"
"        firstPoint[ namesOfVevs[ vevIndexOne ] ] -= firstStepSize\n"
"        firstDifference = ( VPD."
<< PotentialMinimizer::functionFromDictionary << "( inputFunction,\n"
"                                                        vevValues )\n"
"                            - VPD."
<< PotentialMinimizer::functionFromDictionary << "( inputFunction,\n"
"                                                          firstPoint ) )\n"
"        for vevIndexTwo in range( vevIndexOne, numberOfVevs ):\n"
"            secondPoint = vevValues.copy()\n"
"            secondStepSize = stepSize\n"
"            secondPoint[ namesOfVevs[ vevIndexTwo ] ] += secondStepSize\n"
"            doubleOffset = secondPoint.copy()\n"
"            doubleOffset[ namesOfVevs[ vevIndexOne ] ] -= firstStepSize\n"
"            returnHessian[ vevIndexOne ][ vevIndexTwo\n"
"                         ] = ( ( VPD."
<<           PotentialMinimizer::functionFromDictionary << "( inputFunction,\n"
"                                                             secondPoint )\n"
"                               - VPD."
<<           PotentialMinimizer::functionFromDictionary << "( inputFunction,\n"
"                                                             doubleOffset )\n"
"                                 - firstDifference )\n"
"                               / ( firstStepSize * secondStepSize ) )\n"
"            returnHessian[ vevIndexTwo ][ vevIndexOne\n"
"                         ] = returnHessian[ vevIndexOne ][ vevIndexTwo ]\n"
"    return returnHessian\n"
"\n"
<< pointsToTry << " = VTE." << pointsToTry << "\n"
"if ( 0 >= len( " + pointsToTry + " ) ):\n"
"    warningMessage = \"no tree-level extrema found!\"\n"
"    warningMessages.append( warningMessage )\n"
"    print( warningMessage )\n"
"\n"
"foundMinima = []\n"
"foundSaddles = []\n"
"\n"
"minuitObject = minuit.Minuit( VPD."
<<                         PotentialMinimizer::loopCorrectedPotential << " )\n"
"# There is a possibility of MINUIT trying to roll off to infinity but\n"
"# generating an overflow error in calculating the potential before throwing\n"
"# an exception. Here limits are set on the VEVs of one thousand times\n"
"# the energy scale, which is probably way beyond where the potential should\n"
"# be trusted anyway.\n"
"for vevVariable in minuitObject.values:\n"
"    minuitObject.limits[ vevVariable ] = ( -1.0E+3, 1.0E+3 )\n"
"\n"
"\n"
"def VevsHaveCorrectSigns( vevValues ):\n"
"#    print( \"checking signs of \" + VPD."
<<             PotentialMinimizer::userVevsAsMathematica << "( vevValues ) )\n"
"    returnBool = True\n"
"    for positiveVev in VPD." << PotentialMinimizer::vevsTakenPositive << ":\n"
"        if ( -stepSize > vevValues[ positiveVev ] ):\n"
"            returnBool = False\n"
"#    print( \"returning \" + str( returnBool ) )\n"
"    return returnBool\n"
"\n"
"def SteepestDescent( vevValues ):\n"
"    eigensystemOfHessian = numpy.linalg.eigh( NumericalHessian(\n"
"                                                VPD."
<<                          PotentialMinimizer::loopCorrectedPotential << ",\n"
"                                                              vevValues ) )\n"
"    mostNegativeEigenvalueValue = 0.0\n"
"    mostNegativeEigenvalueIndex = 0\n"
"    for eigenvalueIndex in range( len( eigensystemOfHessian[ 0 ] ) ):\n"
"        if ( mostNegativeEigenvalueValue > eigensystemOfHessian[ 0 ][\n"
"                                                       eigenvalueIndex ] ):\n"
"            mostNegativeEigenvalueValue = eigensystemOfHessian[ 0 ][\n"
"                                                          eigenvalueIndex ]\n"
"            mostNegativeEigenvalueIndex = eigenvalueIndex\n"
"    return [ mostNegativeEigenvalueValue,\n"
"             eigensystemOfHessian[ 1 ][ :, mostNegativeEigenvalueIndex ] ]\n"
"\n"
"def TryToMinimize( vevValues ):\n"
"#    print( \"trying to minimize\" + VPD."
<<             PotentialMinimizer::userVevsAsMathematica << "( vevValues ) )\n"
"    global minuitObject\n"
"    global foundMinima\n"
"    global foundSaddles\n"
"    try:\n"
"        minuitObject.values = vevValues.copy()\n"
"        minuitObject.migrad()\n"
"        foundExtremum = { \'potentialDepth\': minuitObject.fval,\n"
"                          \'depthError\': abs( minuitObject.edm ),\n"
"                          \'vevValues\': minuitObject.values.copy() }\n"
"        candidateDescent = SteepestDescent( minuitObject.values )\n"
"        if ( 0.0 > candidateDescent[ 0 ] ):\n"
"            foundSaddles.append( [ foundExtremum,\n"
"                                   list( candidateDescent[ 1 ] ) ] )\n"
"        else:\n"
"            foundMinima.append( foundExtremum )\n"
"    except minuit.MinuitError as minuitError:\n"
"        warningMessage = ( \"PyMinuit had problems starting at \"\n"
"                           + VPD."
<<               PotentialMinimizer::userVevsAsMathematica << "( vevValues )\n"
"                           + \" [minuit.MinuitError: \"\n"
"                           + str( minuitError )\n"
"                           + \"]. PyMinuit stopped at \"\n"
"                          + VPD." << PotentialMinimizer::userVevsAsMathematica
<<                                                  "( minuitObject.values )\n"
"                           + \" with relative depth \"\n"
"               + str( ( VPD." << PotentialMinimizer::energyScaleFourth << "\n"
"                        * VPD."
<<                       PotentialMinimizer::functionFromDictionary << "( VPD."
<<                          PotentialMinimizer::loopCorrectedPotential << ",\n"
"                                                    minuitObject.values ) )\n"
"                      - potentialAtVevOrigin )\n"
"                            + \" at one-loop level and \"\n"
"                 + str( VPD." << PotentialMinimizer::energyScaleFourth << "\n"
"                        * VPD."
<<                       PotentialMinimizer::functionFromDictionary << "( VPD."
<<                              PotentialMinimizer::treeLevelPotential << ",\n"
"                                                    minuitObject.values ) )\n"
"                            + \" at tree level.\"\n"
"               + \" Minuit's estimate of how much deeper it should go is \"\n"
"                 + str( VPD." << PotentialMinimizer::energyScaleFourth << "\n"
"                        * minuitObject.edm ) )\n"
"        warningMessages.append( warningMessage )\n"
"        print( warningMessage )\n"
"\n"
"def DisplacePoint( pointDictionary, displacementList, scaleFactor ):\n"
"    returnDictionary = pointDictionary.copy()\n"
"    for vevIndex in range( numberOfVevs ):\n"
"        returnDictionary[ namesOfVevs[ vevIndex ] ] += ( scaleFactor\n"
"                                           * displacementList[ vevIndex ] )\n"
"    return returnDictionary\n"
"\n"
"for vevValueSet in " << pointsToTry << ":\n"
"    if ( VevsHaveCorrectSigns( vevValueSet ) ):\n"
"        TryToMinimize( vevValueSet )\n"
"\n"
"for saddleSplit in range( VPD." << maximumSaddleSplitsVariableName << " ):\n"
"    if ( 0 < len( foundSaddles ) ):\n"
"        print( \"PyMinuit had to be nudged off \"\n"
"               + str( len( foundSaddles ) ) + \" saddle point(s).\" )\n"
"        " << pointsToTry << " = []\n"
"        for saddlePoint in foundSaddles:\n"
"            " << pointsToTry << ".append( DisplacePoint(\n"
"                                         saddlePoint[ 0 ][ \'vevValues\' ],\n"
"                                                       saddlePoint[ 1 ], \n"
"                                                       stepSize ) )\n"
"            " << pointsToTry << ".append( DisplacePoint(\n"
"                                         saddlePoint[ 0 ][ \'vevValues\' ],\n"
"                                                       saddlePoint[ 1 ], \n"
"                                                       -stepSize ) )\n"
"        foundSaddles = []\n"
"        for vevValueSet in " << pointsToTry << ":\n"
"            if ( VevsHaveCorrectSigns( vevValueSet ) ):\n"
"                TryToMinimize( vevValueSet )\n"
"\n"
"globalMinimumCandidates = foundMinima\n"
"\n"
"# This bit of code adds in any remaining saddle points to\n"
"# globalMinimumCandidates as well as setting up warning messages:\n"
"if ( 0 != len( foundSaddles ) ):\n"
"    warningMessage = ( str( len( foundSaddles ) )\n"
"                       + \" extremum/a with at least one descending or\"\n"
"                     + \" flat direction remained after all nudging:\" )\n"
"    for saddlePoint in foundSaddles:\n"
"        warningMessage += \"\\n \" + VPD."
<<                           PotentialMinimizer::userVevsAsMathematica << "(\n"
"                                        saddlePoint[ 0 ][ \'vevValues\' ] )\n"
"        globalMinimumCandidates.append( saddlePoint[ 0 ] )\n"
"    warningMessages.append( warningMessage )\n"
"    print( warningMessage )\n"
"\n"
"# If the input vacuum is the global minimum, actionValue is set to -1.0.\n"
"# Non-negative values of actionValue indicate the current upper bound on\n"
"# the action after the last approximation. It starts stupidly high\n"
"# (the current age of the Universe corresponds to an action of about 400)\n"
"# so that the first requested bounding estimate will be calculated.\n"
"actionNeedsToBeCalculated = False\n"
"stabilityVerdict = \"error\"\n"
"actionValue = 1000000.0\n"
"actionType = \"error\"\n"
"tunnelingTime = -2.0\n"
"\n"
"if ( 0 >= len( globalMinimumCandidates ) ):\n"
"    warningMessage = \"no 1-loop-level extrema found!\"\n"
"    warningMessages.append( warningMessage )\n"
"    print( warningMessage )\n"
"\n"
"# The result is assumed metastable unless found otherwise.\n"
"stabilityVerdict = \"metastable\"\n"
"givenInputAsArray = VPD." << PotentialMinimizer::vevDictionaryToArray
<<                     "( VPD." << PotentialMinimizer::inputVevsPoint << " )\n"
"minuitObject.values = VPD." << PotentialMinimizer::inputVevsPoint
<<                                                                  ".copy()\n"
"foundSaddles = []\n"
"try:\n"
"    minuitObject.migrad()\n"
"except minuit.MinuitError as minuitError:\n"
"    warningMessage = ( \"PyMinuit had problems starting at input VEVs! \"\n"
"                     + VPD"
<<               PotentialMinimizer::userVevsAsMathematica << "( vevValues )\n"
"                       + \" [minuit.MinuitError: \"\n"
"                       + str( minuitError )\n"
"                       + \"]. PyMinuit stopped at \"\n"
"                       + VPD." << PotentialMinimizer::userVevsAsMathematica
<<                                                  "( minuitObject.values )\n"
"                       + \" with relative depth \"\n"
"               + str( ( VPD." << PotentialMinimizer::energyScaleFourth << "\n"
"                        * VPD."
<<                       PotentialMinimizer::functionFromDictionary << "( VPD."
<<                          PotentialMinimizer::loopCorrectedPotential << ",\n"
"                                                    minuitObject.values ) )\n"
"                      - potentialAtVevOrigin )\n"
"                            + \" at one-loop level and \"\n"
"                 + str( VPD." << PotentialMinimizer::energyScaleFourth << "\n"
"                        * VPD."
<<                       PotentialMinimizer::functionFromDictionary << "( VPD."
<<                              PotentialMinimizer::treeLevelPotential << ",\n"
"                                                    minuitObject.values ) )\n"
"                        + \" at tree level.\"\n"
"               + \" Minuit's estimate of how much deeper it should go is \"\n"
"                 + str( VPD." << PotentialMinimizer::energyScaleFourth << "\n"
"                        * minuitObject.edm ) )\n"
"    warningMessages.append( warningMessage )\n"
"    print( warningMessage )\n"
"if ( 0 != len( foundSaddles ) ):\n"
"    warningMessage = \"PyMinuit rolled from input VEVs to a saddle point!\"\n"
"    warningMessages.append( warningMessage )\n"
"    print( warningMessage )\n"
"rolledInputAsDictionary = minuitObject.values.copy()\n"
"# Occasionally degenerate minima with different signs, that are equivalent\n"
"# under a gauge transformation, appear, & by the nature of the algorithm,\n"
"# MINUIT might roll closer to the sign-flipped minimum than how close it\n"
"# rolls to the input minimum. In such cases, the (extremely long) tunneling\n"
"# time to the gauge-equivalent minimum would be calculated. To avoid this,\n"
"# the input minimum is conservatively taken to have a depth equal to that\n"
"# found by MINUIT minus twice MINUIT's error.\n"
"rolledInputDepth = minuitObject.fval - abs( 2.0 * minuitObject.edm )\n"
"# The input VEV point is taken as the global minimum to begin with:\n"
"globalMinimumPointAsDictionary = rolledInputAsDictionary.copy()\n"
"globalMinimumDepthError = abs( minuitObject.edm )\n"
"globalMinimumDepthValue = rolledInputDepth\n"
"rolledInputAsArray = VPD." << PotentialMinimizer::vevDictionaryToArray
<<                                              "( rolledInputAsDictionary )\n"
"rollDistanceSquared = numpy.sum( ( rolledInputAsArray\n"
"                                   - givenInputAsArray )**2 )\n"
"if ( ( rollingToleranceSquared * numpy.sum( givenInputAsArray**2 ) )\n"
"     < rollDistanceSquared ):\n"
"    warningMessage = (\n"
"                  \"PyMinuit rolled quite far from the input VEVs! (from \"\n"
"             + str( VPD." << PotentialMinimizer::userVevsAsMathematica
<<        "( VPD." << PotentialMinimizer::inputVevsPoint << " ) ) + \" to \"\n"
"             + str( VPD." << PotentialMinimizer::userVevsAsMathematica
<<                                      "( minuitObject.values ) ) + \")\" )\n"
"    warningMessages.append( warningMessage )\n"
"    print( warningMessage )\n"
"rolledInputLengthSquared = numpy.sum( rolledInputAsArray**2 )\n"
"\n"
"outputFile = open( \"./Vevacious_loop-corrected_results.txt\", \"w\" )\n"
"\n"
"for globalMinimumCandidate in globalMinimumCandidates:\n"
"    outputFile.write( str( globalMinimumCandidate[ \'potentialDepth\' ]\n"
"                           * VPD."
<<                     PotentialMinimizer::energyScaleFourth << " ) + \", \"\n"
"                      + VPD." << PotentialMinimizer::userVevsAsMathematica
<<                                                                        "(\n"
"                      globalMinimumCandidate[ \'vevValues\' ] ) + \"\\n\" )\n"
"    if ( globalMinimumDepthValue\n"
"         > globalMinimumCandidate[ \'potentialDepth\' ] ):\n"
"        actionNeedsToBeCalculated = True\n"
"        globalMinimumPointAsDictionary = globalMinimumCandidate[\n"
"                                                     \'vevValues\' ].copy()\n"
"        globalMinimumDepthValue = globalMinimumCandidate[\n"
"                                                       \'potentialDepth\' ]\n"
"        globalMinimumDepthError = globalMinimumCandidate[ \'depthError\' ]\n"
"\n"
"outputFile.close()\n"
"\n"
"numericallyDegenerateCandidates = []\n"
"for globalMinimumCandidate in globalMinimumCandidates:\n"
"    if ( ( globalMinimumCandidate[ \'depthError\' ]\n"
"           + globalMinimumDepthError )\n"
"         > ( globalMinimumCandidate[ \'potentialDepth\' ]\n"
"             - globalMinimumDepthValue ) ):\n"
"       numericallyDegenerateCandidates.append(\n"
"                                            globalMinimumCandidate.copy() )\n"
"\n"
"shortestDistanceSquared = -1.0\n"
"for closeEnoughMinimum in numericallyDegenerateCandidates:\n"
"    if ( rolledInputDepth > closeEnoughMinimum[ \'potentialDepth\' ] ):\n"
"        closeEnoughMinimumAsArray = VPD."
<<                            PotentialMinimizer::vevDictionaryToArray << "(\n"
"                                      closeEnoughMinimum[ \'vevValues\' ] )\n"
"        distanceSquaredToInput = numpy.sum( ( closeEnoughMinimumAsArray\n"
"                                              - rolledInputAsArray )**2 )\n"
"# shortestDistanceSquared is negative before the 1st\n"
"# candidate minimum has been checked for its distance from the input\n"
"# VEVs, so this 1st \"if\" ensures that the 1st candidate minimum\n"
"# is taken as the best minimum so far automatically.\n"
"        if ( shortestDistanceSquared < 0.0 ):\n"
"            shortestDistanceSquared = ( distanceSquaredToInput + 1.0 )\n"
"        if ( distanceSquaredToInput < shortestDistanceSquared ):\n"
"            shortestDistanceSquared = distanceSquaredToInput\n"
"            globalMinimumPointAsDictionary = closeEnoughMinimum[\n"
"                                                     \'vevValues\' ].copy()\n"
"\n"
"globalMinimumPointAsArray = VPD." << PotentialMinimizer::vevDictionaryToArray
<<                                                                        "(\n"
"                                           globalMinimumPointAsDictionary )\n"
"\n"
"# if the global minimum is sufficiently far away from the input VEVs that\n"
"# it is unlikely that is was just that MINUIT did not roll exactly to where\n"
"# there should have been a minimum at the input VEVs...\n"
"distanceSquaredFromInputToGlobalMinimum = numpy.sum( (\n"
"                      globalMinimumPointAsArray - rolledInputAsArray )**2 )\n"
"\n"
"# If the input vacuum is the global minimum, actionValue is set to -1.0.\n"
"# Non-negative values of actionValue indicate the current upper bound on\n"
"# the action after the last approximation. It starts stupidly high\n"
"# (the current age of the Universe corresponds to an action of about 400)\n"
"# so that the first requested bounding estimate will be calculated.\n"
"if ( ( rollingToleranceSquared * rolledInputLengthSquared )\n"
"     > distanceSquaredFromInputToGlobalMinimum ):\n"
"    stabilityVerdict = \"stable\"\n"
"    actionValue = -1.0\n"
"    tunnelingTime = -1.0\n"
"    actionType = \"unnecessary\"\n"
"    actionNeedsToBeCalculated = False\n"
"\n"
"if ( actionNeedsToBeCalculated\n"
"     and ( 0.0 < VPD." << quarticActionBoundVariableName << " )\n"
"     and ( actionValue > VPD." << quarticActionBoundVariableName << " ) ):\n"
"# These numbers come from fitting a quartic polynomial to have its minima\n"
"# at the (rolled) input VEVs and the global minimum, also having the\n"
"# correct value half-way between them. then the cubic term is thrown away,\n"
"# so the quartic now sits atop the full potential, snugly at the input\n"
"# minimum, and the coefficients matched to those given by Sidney Coleman in\n"
"# PRD vol.15, num.10, 15th May 1977 (i.e. {\\mu}^{2}, here muSquared, and\n"
"# \\lambda, here lambdaValue. Then Coleman\'s epsilon, here epsilonValue,\n"
"# is identified as the difference in depths at the minima. This should give\n"
"# a very conservative upper bound on the tunneling time.\n"
"    halfwayPoint = ( ( globalMinimumPointAsArray\n"
"                       + rolledInputAsArray ) * 0.5 )\n"
"    halfwayDepth = VPD." << PotentialMinimizer::functionFromArray << "( VPD."
<<                          PotentialMinimizer::loopCorrectedPotential << ",\n"
"                                          halfwayPoint )\n"
"    muSquared = ( ( 2.0 * ( 16.0 * halfwayDepth\n"
"                            - 11.0 * rolledInputDepth\n"
"                            - 5.0 * globalMinimumDepthValue ) )\n"
"                  / distanceSquaredFromInputToGlobalMinimum )\n"
"    lambdaValue = ( ( 8.0 * ( 2.0 * halfwayDepth - rolledInputDepth\n"
"                              - globalMinimumDepthValue ) )\n"
"                    / distanceSquaredFromInputToGlobalMinimum**2 )\n"
"    epsilonValue = ( rolledInputDepth - globalMinimumDepthValue )\n"
"    actionValue = ( ( math.pi * math.pi * muSquared**6 )\n"
"                    / ( 6.0 * lambdaValue**4 * epsilonValue**3 ) )\n"
"    actionType = \"quartic_bound\"\n"
"    if ( actionValue < VPD." << quarticActionBoundVariableName << " ):\n"
"        stabilityVerdict = \"unstable\"\n"
"        actionNeedsToBeCalculated = False\n"
"\n"
"# The resolution of the tunneling path needs to be set\n"
"# (low-ish by default for speed):\n"
"tunnelingResolution = 20\n"
"\n"
"if ( actionNeedsToBeCalculated\n"
"     and ( ( 0.0 < VPD." << directActionBoundVariableName << " )\n"
"           or ( 0.0 < VPD." << deformedActionBoundVariableName << " ) ) ):\n"
"    firstStepPoint = ( ( rolledInputAsArray\n"
"                         * ( 1.0 - ( 1.0 / tunnelingResolution ) ) )\n"
"                       + ( globalMinimumPointAsArray\n"
"                           * ( 1.0 / tunnelingResolution ) ) )\n"
"    firstStepDepth = VPD." << PotentialMinimizer::functionFromArray
<<              "( VPD." << PotentialMinimizer::loopCorrectedPotential << ",\n"
"                                            firstStepPoint )\n"
"    if ( rolledInputDepth >= firstStepDepth ):\n"
"        actionType = \"barrier_smaller_than_resolution\"\n"
"        actionValue = 0.0\n"
"        tunnelingTime = 0.0\n"
"        stabilityVerdict = \"unstable\"\n"
"        actionNeedsToBeCalculated = False\n"
"        warningMessage = ( \"Energy barrier from input VEVs to global\"\n"
"                           + \" minimum less than resolution of tunneling\"\n"
"                           + \" path!\" )\n"
"        warningMessages.append( warningMessage )\n"
"        print( warningMessage )\n"
"\n"
"if ( actionNeedsToBeCalculated\n"
"     and ( ( 0.0 < VPD." << directActionBoundVariableName << " )\n"
"           or ( 0.0 < VPD." << deformedActionBoundVariableName << " ) ) ):\n"
"    import sys\n"
"    sys.path.append( VPD." << pathToCosmotransitionsVariableName << " )\n"
"    import pathDeformation as CPD\n"
"\n"
"    arrayOfArrays = numpy.array( [ globalMinimumPointAsArray.copy(),\n"
"                                   rolledInputAsArray.copy() ] )\n"
"\n"
"    def PotentialFromArray( pointAsArray ):\n"
"        return VPD." << PotentialMinimizer::functionFromArray << "( VPD."
<<                          PotentialMinimizer::loopCorrectedPotential << ",\n"
"                                                             pointAsArray )\n"
"\n"
"    def PotentialFromMatrix( arrayOfArrays ):\n"
"        if ( ( numberOfVevs, ) == arrayOfArrays.shape ):\n"
"            return PotentialFromArray( arrayOfArrays )\n"
"        elif ( ( len( arrayOfArrays ), numberOfVevs )\n"
"               == arrayOfArrays.shape ):\n"
"            returnArray = numpy.zeros( len( arrayOfArrays ) )\n"
"            for whichIndex in range( len( arrayOfArrays ) ):\n"
"                returnArray[ whichIndex ] = PotentialFromArray(\n"
"                                              arrayOfArrays[ whichIndex ] )\n"
"            return returnArray\n"
"        else:\n"
"            return None\n"
"\n"
"    def GradientFromArray( pointAsArray ):\n"
"        potentialAtPoint = PotentialFromArray( pointAsArray )\n"
"        gradientArray = numpy.zeros( len( pointAsArray ) )\n"
"        for whichField in range( len( pointAsArray ) ):\n"
"            displacedPoint = pointAsArray.copy()\n"
"            displacedPoint[ whichField ] += stepSize\n"
"            gradientArray[ whichField ] = (\n"
"                                     ( PotentialFromArray( displacedPoint )\n"
"                                       - potentialAtPoint ) / stepSize )\n"
"        return gradientArray\n"
"\n"
"    def GradientFromMatrix( arrayOfArrays ):\n"
"        if ( ( numberOfVevs, ) == arrayOfArrays.shape ):\n"
"            return GradientFromArray( arrayOfArrays )\n"
"        elif ( ( len( arrayOfArrays ), numberOfVevs )\n"
"               == arrayOfArrays.shape ):\n"
"            returnMatrix = arrayOfArrays.copy()\n"
"            for whichIndex in range( len( arrayOfArrays ) ):\n"
"                returnMatrix[ whichIndex ] = GradientFromArray(\n"
"                                              arrayOfArrays[ whichIndex ] )\n"
"            return returnMatrix\n"
"        else:\n"
"            return None\n"
"\n"
"    if ( ( 0.0 < VPD." << directActionBoundVariableName << " )\n"
"         and ( actionValue > VPD." << directActionBoundVariableName
<<                                                                    " ) ):\n"
"        quickTunneler = CPD.fullTunneling( V = PotentialFromMatrix,\n"
"                                           dV = GradientFromMatrix,\n"
"                                           phi = arrayOfArrays,\n"
"                                           quickTunneling = True,\n"
"                                           npoints = tunnelingResolution )\n"
"        quickTunneler.doQuickTunnel()\n"
"        quickTunneler.tunnel1D( xtol = 1e-4, phitol = 1e-6 )\n"
"        actionValue = quickTunneler.findAction()\n"
"        actionType = \"direct_path_bound\"\n"
"        if( actionValue < VPD." << directActionBoundVariableName << " ):\n"
"            stabilityVerdict = \"unstable\"\n"
"            actionNeedsToBeCalculated = False\n"
"\n"
"    if ( actionNeedsToBeCalculated\n"
"         and ( 0.0 < VPD." << deformedActionBoundVariableName << " )\n"
"         and ( actionValue > VPD." + deformedActionBoundVariableName
<<                                                                    " ) ):\n"
"        fullTunneler = CPD.fullTunneling( V = PotentialFromMatrix,\n"
"                                          dV = GradientFromMatrix,\n"
"                                          phi = arrayOfArrays,\n"
"                                          quickTunneling = True,\n"
"                                          npoints = tunnelingResolution )\n"
"        fullTunneler.run()\n"
"        actionValue = fullTunneler.findAction()\n"
"        actionType = \"full_deformed_path\"\n"
"        if ( actionValue < VPD." << deformedActionBoundVariableName << " ):\n"
"            stabilityVerdict = \"unstable\"\n"
"            actionNeedsToBeCalculated = False\n"
"\n"
"# No matter if there were serious errors or not, an output file is written:\n"
"outputFile = open( VPD." << resultsFilenameVariableName << ", \"w\" )\n"
"outputFile.write( \"<Vevacious_result>\\n  <stability> \"\n"
"                      + stabilityVerdict\n"
"              + \" </stability>\\n  <global_minimum   relative_depth=\\\"\"\n"
"                      + str( ( globalMinimumDepthValue * VPD."
<<                                        PotentialMinimizer::energyScaleFourth
<<                                  " ) - potentialAtVevOrigin ) + \"\\\" \"\n"
"                      + VPD." << PotentialMinimizer::userVevsAsXml
<<                                       "( globalMinimumPointAsDictionary )\n"
"                      + \" />\\n  <input_minimum   relative_depth=\\\"\"\n"
"                      + str( ( rolledInputDepth * VPD."
<<                                        PotentialMinimizer::energyScaleFourth
<<                                  " ) - potentialAtVevOrigin ) + \"\\\" \"\n"
"                      + VPD." << PotentialMinimizer::userVevsAsXml
<<                                              "( rolledInputAsDictionary )\n"
"                      + \" />\" )\n"
"if ( 0.0 <= actionValue ):\n"
"# We don't want an overflow error when calculating the tunneling time.\n"
"    if ( 1000.0 < actionValue ):\n"
"        actionValue = 1000.0\n"
"    tunnelingTime = ( math.exp( 0.25 * actionValue ) * 1.0e-43 )\n"
"    if ( 1000000.0 < tunnelingTime ):\n"
"        tunnelingTime = 1000000.0\n"
"# The tunneling time is capped at a million times the current age of the\n"
"# known Universe.\n"
"# This code assumes that the age of the Universe is 10^(44)/TeV, and that\n"
"# the solitonic solution factor (Sidney Coleman\'s \"A\") is (0.1 TeV)^4.\n"
"outputFile.write( \"\\n  <lifetime  action_calculation=\\\"\" + actionType\n"
"                  + \"\\\" > \" + str( tunnelingTime ) + \" </lifetime>\" )\n"
"# Each warning is printed as an XML element:\n"
"for warningMessage in warningMessages:\n"
"    outputFile.write( \"\\n  <warning>\\n  \" + warningMessage\n"
"                      + \"\\n  </warning>\" )\n"
"outputFile.write( \"\\n</Vevacious_result>\\n\" )\n"
"outputFile.close()\n";
    outputFile.close();
  }

} /* namespace Vevacious */
