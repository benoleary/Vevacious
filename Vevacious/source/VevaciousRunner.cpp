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
  std::string const VevaciousRunner::vevaciousVersion( "1.1.01" );
  std::string const VevaciousRunner::vevaciousDocumentation(
                                 "arXiv:1307.1477, arXiv:1405.7376 (hep-ph)" );
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
    resultsFilename( "MyResult.vout" ),
    imaginaryTolerance( 0.000001 ),
    saddleNudgeList( 2,
                     1.0 ),
    rollingTolerance( 0.2 ),
    pathToCosmotransitions( "./" ),
    lifetimeThreshold( 0.1 ),
    survivalThreshold( 0.1 ),
    shouldTunnel( "True" ),
    tunnelThermally( "True" ),
    deformTunnelPaths( "True" )
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
                                             "./MyResult.vout" ) ),
    imaginaryTolerance( BOL::StringParser::stringToDouble(
                                 argumentParser.fromTag( "imaginary_tolerance",
                                                         "0.0000001" ) ) ),
    saddleNudgeList( 2,
                     1.0 ),
    rollingTolerance( BOL::StringParser::stringToDouble(
                                      argumentParser.fromTag( "roll_tolerance",
                                                              "0.1" ) ) ),
    pathToCosmotransitions( argumentParser.fromTag( "ct_path",
                                                    "./CosmoTransitions/" ) ),
    lifetimeThreshold( BOL::StringParser::stringToDouble(
                                  argumentParser.fromTag( "lifetime_threshold",
                                                               "0.1" ) ) ),
    survivalThreshold( BOL::StringParser::stringToDouble(
                          argumentParser.fromTag( "thermal_survival_threshold",
                                                          "0.1" ) ) ),
    shouldTunnel( argumentParser.fromTag( "should_tunnel",
                                          "True" ) ),
    tunnelThermally( argumentParser.fromTag( "tunnel_thermally",
                                             "True" ) ),
    deformTunnelPaths( argumentParser.fromTag( "deform_tunnel_paths",
                                               "True" ) )
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "shouldTunnel = \"" << shouldTunnel << "\""
    << std::endl
    << "tunnelThermally = \"" << tunnelThermally << "\""
    << std::endl
    << "deformTunnelPaths = \"" << deformTunnelPaths << "\"";
    std::cout << std::endl;/**/

    setMinuitNudgesOffSaddlePoints( argumentParser.fromTag( "saddle_nudges",
                                                            "1.0, 1.0" ) );
    std::string
    inputMaximumSaddleNudges( argumentParser.fromTag( "max_saddle_nudges",
                                                      "" ) );
    if( !(inputMaximumSaddleNudges.empty()) )
    {
      setMaximumMinuitNudgesOffSaddlePoints( inputMaximumSaddleNudges );
    }
    BOL::AsciiXmlParser xmlParser;
    std::string modelFilename( argumentParser.fromTag( "model_file",
                                                       "./MyModel.vin" ) );
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
"import scipy\n"
"import sys\n"
"import minuit\n"
"import time\n"
"\n"
"# This file was automatically generated by Vevacious\n"
"# version " << vevaciousVersion << ".\n"
"# It contains 2 classes & several global variables & functions.\n"
"# These appear in the order:\n"
"# global variables and functions\n"
"# FunctionScaler class\n"
"# Vevacious class\n"
"\n"
"\n"
"\n"
"outputFile = \"" << resultsFilename << "\"\n"
"lnOfThermalIntegrationFactor = 244.53\n"
"energyScale = " << sarahInterpreter.getSlhaScale() << "\n"
"energyScaleFourth = ( energyScale**4 )\n"
"inverseScale = ( 1.0 / energyScale )\n"
"inverseScaleSquared = ( 1.0 / ( energyScale * energyScale ) )\n"
"inverseScaleFourthed = ( 1.0 / energyScaleFourth )\n"
"numericalStepSize = 1.0\n"
"exponentCutOff = math.log( 0.5 * sys.float_info.max )\n"
"\n"
"\n"
"# First some general functions:\n"
"\n"
"pathToCosmotransitions = \"" << pathToCosmotransitions << "\"\n"
"sys.path.append( pathToCosmotransitions )\n"
"import pathDeformation as CTPD\n"
"ctVersionString = getattr( CTPD, \"__version__\", \"1\" )\n"
"ctMajorVersion = int( ctVersionString.split( \'.\' )[ 0 ] )\n"
"import finiteT as CTFT\n"
"\n# Unfortunately due to UTF32 issues, PyMinuit is restricted to\n"
"# single-character variable names, so Vevacious internally renames the\n"
"# fields:\n"
"# " << sarahInterpreter.getHumanReadableVevNameMap() << "\n"
"namesOfFields = [ "
<<             sarahInterpreter.getInternalVevNamesAsQuotedCharList() << " ]\n"
"numberOfFields = len( namesOfFields )\n"
"internalFieldNamesToUserNames = { ";
    std::vector< char > const&
    vevsOrderedBySolution( sarahInterpreter.getVevsOrderedBySolutions() );
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        outputFile << ",\n                                  ";
      }
      outputFile << "\'" << *whichName << "\': \""
      << sarahInterpreter.getUserVevName( *whichName ) << "\"";
    }
    outputFile << " }\n"
"vevsTakenPositive = " << sarahInterpreter.getPositiveInternalVevs() << "\n"
"fieldOrigin = { ";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        outputFile << ", ";
      }
      outputFile << "\'" << *whichName << "\': 0.0";
    }
    outputFile << " }\n"
"originArray = [ ";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        outputFile << ", ";
      }
      outputFile << "0.0";
    }
    outputFile << " ]\n"
"dsbInput = { ";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        outputFile << ",\n             ";
      }
      outputFile << "\'" << *whichName << "\': ( "
      << sarahInterpreter.getInputVevValue( *whichName ) << " )";
    }
    outputFile << " }\n"
"\n"
"\n"
"def FieldDictionaryToArray( fieldValueDictionary ):\n"
"    return numpy.array( [ ";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        outputFile << ",\n                          ";
      }
      outputFile << "fieldValueDictionary[ \'" << *whichName << "\' ]";
    }
    outputFile << " ] )\n"
"\n"
"\n"
"def FieldArrayToDictionary( fieldValueArray ):\n"
"    return { ";
    for( unsigned int whichVev( 0 );
         vevsOrderedBySolution.size() > whichVev;
         ++whichVev )
    {
      if( 0 < whichVev )
      {
        outputFile << ",\n             ";
      }
      outputFile
      << "\'" << vevsOrderedBySolution[ whichVev ] << "\': fieldValueArray[ "
      << whichVev << " ]";
    }
    outputFile << " }\n"
"\n"
"\n"
"def UserFieldsAsMathematica( fieldValueDictionary ):\n"
"    return ( \"{ ";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        outputFile << ", ";
      }
      outputFile << sarahInterpreter.getUserVevName( *whichName )
      << " -> ( \"\n             + str( fieldValueDictionary[ \'"
      << *whichName << "\' ] )\n                 + \" )";
    }
    outputFile << " }\" )\n"
"\n"
"\n"
"def UserFieldsAsXml( fieldValueDictionary ):\n"
"    return ( \"";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        outputFile << " ";
      }
      outputFile << sarahInterpreter.getUserVevName( *whichName )
      << "=\\\"\"\n               + str( fieldValueDictionary[ \'"
      << *whichName << "\' ] ) + \"\\\"";
    }
    outputFile << "\" )\n"
"\n"
"\n"
"def FunctionFromDictionary( FunctionFromArguments,\n"
"                            fieldValueDictionary,\n"
"                            givenTemperatureValue = 0.0 ):\n"
"    return FunctionFromArguments( ";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      outputFile << *whichName << " = fieldValueDictionary[ \'" << *whichName
      << "\' ], \n                                  ";
    }
    outputFile << "temperatureValue = givenTemperatureValue )\n"
"\n"
"\n"
"def FunctionFromArray( FunctionFromArguments,\n"
"                       fieldArray,\n"
"                       givenTemperatureValue = 0.0 ):\n"
"    return FunctionFromArguments( ";
    for( unsigned int whichVev( 0 );
         vevsOrderedBySolution.size() > whichVev;
         ++whichVev )
    {
      outputFile << vevsOrderedBySolution[ whichVev ] << " = fieldArray[ "
      << whichVev << " ],\n                                  ";
    }
    outputFile << "temperatureValue = givenTemperatureValue )\n"
"\n"
"\n"
"def LengthSquaredOfMinimum( givenMinimum ):\n"
"    return numpy.sum( (FieldDictionaryToArray( givenMinimum[\n"
"                                                  \"FieldValues\" ] ))**2 )\n"
"\n"
"\n"
"# MINUIT\'s hesse() function assumes that it is already at a\n"
"# minimum, but we need to check whether it actually stopped at a\n"
"# saddle point, so we need to work out the Hessian matrix ourselves.\n"
"def NumericalHessian( inputFunction, fieldValueDictionary ):\n"
"    numberOfFields = len( fieldValueDictionary )\n"
"    returnHessian = [ [ 0.0 for fieldIndexOne in range( numberOfFields ) ]\n"
"                             for fieldIndexTwo in range( numberOfFields ) ]\n"
"    for fieldIndexOne in range( numberOfFields ):\n"
"        firstPoint = fieldValueDictionary.copy()\n"
"        firstPoint[ namesOfFields[ fieldIndexOne ] ] -= numericalStepSize\n"
"        firstDifference = ( FunctionFromDictionary( inputFunction,\n"
"                                                    fieldValueDictionary )\n"
"                            - FunctionFromDictionary( inputFunction,\n"
"                                                      firstPoint ) )\n"
"        for fieldIndexTwo in range( fieldIndexOne, numberOfFields ):\n"
"            secondPoint = fieldValueDictionary.copy()\n"
"            secondPoint[ namesOfFields[ fieldIndexTwo\n"
"                                                   ] ] += numericalStepSize\n"
"            doubleOffset = secondPoint.copy()\n"
"            doubleOffset[ namesOfFields[ fieldIndexOne\n"
"                                                   ] ] -= numericalStepSize\n"
"            returnHessian[ fieldIndexOne ][ fieldIndexTwo\n"
"                             ] = ( ( FunctionFromDictionary( inputFunction,\n"
"                                                             secondPoint )\n"
"                                   - FunctionFromDictionary( inputFunction,\n"
"                                                             doubleOffset )\n"
"                                     - firstDifference )\n"
"                              / ( numericalStepSize * numericalStepSize ) )\n"
"            returnHessian[ fieldIndexTwo ][ fieldIndexOne\n"
"                        ] = returnHessian[ fieldIndexOne ][ fieldIndexTwo ]\n"
"    return returnHessian\n"
"\n"
"\n"
"def MassSquaredCorrections( massSquaredValues,\n"
"                            overallFactor,\n"
"                            subtractionConstant ):\n"
"    summedCorrection = 0.0\n"
"    for massSquaredValue in massSquaredValues:\n"
"        massSquaredMagnitude = abs( massSquaredValue )\n"
"        if( massSquaredMagnitude > 1.0 ):\n"
"            summedCorrection += ( massSquaredMagnitude\n"
"                                  * massSquaredMagnitude\n"
"                                  * ( math.log( massSquaredMagnitude\n"
"                                                * inverseScaleSquared )\n"
"                                      - subtractionConstant ) )\n"
"    return ( overallFactor * summedCorrection )\n"
"\n"
"loopFactor = ( 1.0 / ( 64.0 * math.pi * math.pi ) )\n"
"\n"
"\n"
"# The following are approximations to the J_{+} & J_{-} functions:\n"
"# J_{+}( r ) = integral from 0 to infinity by dx of\n"
"# x^2 ln( 1 - exp( -sqrt( x^2 + r^2 ) ) )\n"
"# J_{-}( r ) = integral from 0 to infinity by dx of\n"
"# x^2 ln( 1 + exp( -sqrt( x^2 + r^2 ) ) )\n"
"# A bosonic degree of freedom (real scalar, so complex scalars count\n"
"# twice!) uses +J_{+}, while a fermionic degree of freedom uses -J_{-} (and\n"
"# counts 2 degrees of freedom per Weyl fermion, so e.g. tops at the SM\n"
"# minimum count 12 times = 3 colors * 2 chiralities * 2 for complex\n"
"# conjugation). The function ThermalFunction returns -J_{-}, so should be\n"
"# added to the potential, not subtracted, i.e. ThermalFunction\'s return\n"
"# value should always be added, multiplied by the positive of the number of\n"
"# degrees of freedom, regardless of whether it is for a boson or a fermion,\n"
"# as it internally accounts for the sign of the spin statistics, but not\n"
"# the absolute number of degrees of freedom!\n"
"\n"
"explicitEulerGamma = 0.577215661901532\n"
"bosonLowA = ( ( -1.0 * (math.pi)**4 ) / 45.0 )\n"
"bosonLowB = ( ( 1.0 * (math.pi)**2 ) / 12.0 )\n"
"bosonLowC = ( ( -1.0 * math.pi ) / 6.0 )\n"
"bosonLowD = ( -1.0 / 32.0 )\n"
"bosonLowLog = ( 1.5 - ( 2.0 * explicitEulerGamma )\n"
"                    + ( 2.0 * math.log( 4.0 * math.pi ) ) )\n"
"\n"
"def LowRatioBosonApproximation( massOverTemperature ):\n"
"    returnValue = bosonLowA\n"
"    if ( massOverTemperature > 0.0 ):\n"
"        returnValue += ( massOverTemperature**2\n"
"                         * ( bosonLowB\n"
"                             + ( massOverTemperature\n"
"                                 * ( bosonLowC\n"
"                                     + ( bosonLowD\n"
"                                         * massOverTemperature\n"
"                                        * ( math.log( massOverTemperature )\n"
"                                             - bosonLowLog ) ) ) ) ) )\n"
"    return returnValue\n"
"\n"
"\n"
"fermionLowA = ( ( -7.0 * (math.pi)**4 ) / 360.0 )\n"
"fermionLowB = ( ( 1.0 * (math.pi)**2 ) / 24.0 )\n"
"fermionLowD = ( 1.0 / 32.0 )\n"
"fermionLowLog = ( 1.5 - ( 2.0 * explicitEulerGamma )\n"
"                      + ( 2.0 * math.log( 1.0 * math.pi ) ) )\n"
"\n"
"def LowRatioFermionApproximation( massOverTemperature ):\n"
"    returnValue = fermionLowA\n"
"    if ( massOverTemperature > 0.0 ):\n"
"        returnValue += ( massOverTemperature**2\n"
"                         * ( fermionLowB\n"
"                         + ( fermionLowD\n"
"                             * massOverTemperature**2\n"
"                             * ( math.log( massOverTemperature )\n"
"                                           - fermionLowLog ) ) ) )\n"
"    return returnValue\n"
"\n"
"\n"
"# In the approximation we work to, the high-ratio (low-temperature)\n"
"# approximation is the same for both bosons and fermions.\n"
"def HighRatioApproximation( massOverTemperature ):\n"
"    return ( -1.0 * massOverTemperature**2\n"
"             * scipy.special.kn( 2, massOverTemperature ) )\n"
"\n"
"# The function for intermediate ratios is chosen as a linear interpolation.\n"
"# The low-ratio approximations are surprisingly accurate up to 1.0!\n"
"lowerKinkRatio = 1.0\n"
"lowerKinkBoson = LowRatioBosonApproximation( lowerKinkRatio )\n"
"lowerKinkFermion = LowRatioFermionApproximation( lowerKinkRatio )\n"
"# The high-ratio approximation is not as good, though is within 3% by 2.5:\n"
"upperKinkRatio = 2.5\n"
"upperKinkValue = HighRatioApproximation( upperKinkRatio )\n"
"slopeBoson = ( ( upperKinkValue - lowerKinkBoson )\n"
"               / ( upperKinkRatio - lowerKinkRatio ) )\n"
"interceptBoson = ( lowerKinkBoson - ( slopeBoson * lowerKinkRatio ) )\n"
"slopeFermion = ( ( upperKinkValue - lowerKinkFermion )\n"
"                 / ( upperKinkRatio - lowerKinkRatio ) )\n"
"interceptFermion = ( lowerKinkFermion - ( slopeFermion * lowerKinkRatio ) )\n"
"\n"
"def IntermediateRatioBosonApproximation( massOverTemperature ):\n"
"    return ( ( slopeBoson * massOverTemperature ) + interceptBoson )\n"
"\n"
"def IntermediateRatioFermionApproximation( massOverTemperature ):\n"
"    return ( ( slopeFermion * massOverTemperature ) + interceptFermion )\n"
"\n"
"thermalInterpolationResolution = 100\n"
"# negativeBosonThermalFunctionValues is a list of length 102 (entries from\n"
"# [0] to [101]) with values for negative m^2/T^2 ratios for\n"
"# J_{+}( m^2/T^2 ) calculated with Mathematica.\n"
"negativeBosonThermalFunctionValues = [ bosonLowA,\n"
"-2.4584,\n"
"-2.70421,\n"
"-2.90927,\n"
"-3.07928,\n"
"-3.21499,\n"
"-3.32186,\n"
"-3.39489,\n"
"-3.45454,\n"
"-3.4678,\n"
"-3.46382,\n"
"-3.44076,\n"
"-3.39106,\n"
"-3.32473,\n"
"-3.23573,\n"
"-3.13883,\n"
"-2.99525,\n"
"-2.84771,\n"
"-2.68783,\n"
"-2.5113,\n"
"-2.31497,\n"
"-2.10299,\n"
"-1.86926,\n"
"-1.62777,\n"
"-1.37701,\n"
"-1.13041,\n"
"-0.84233,\n"
"-0.564025,\n"
"-0.251039,\n"
"0.0627693,\n"
"0.394232,\n"
"0.72439,\n"
"1.07055,\n"
"1.40908,\n"
"1.74135,\n"
"2.13272,\n"
"2.48344,\n"
"2.88792,\n"
"3.25922,\n"
"3.64481,\n"
"4.04158,\n"
"4.45061,\n"
"4.82576,\n"
"5.24276,\n"
"5.6749,\n"
"6.09947,\n"
"6.51899,\n"
"6.92073,\n"
"7.21596,\n"
"7.76782,\n"
"8.17944,\n"
"8.62316,\n"
"9.05358,\n"
"9.44169,\n"
"9.85741,\n"
"10.2806,\n"
"10.6692,\n"
"11.0742,\n"
"11.2484,\n"
"11.7149,\n"
"12.3143,\n"
"12.6849,\n"
"13.0811,\n"
"13.3769,\n"
"13.7479,\n"
"14.2238,\n"
"14.4005,\n"
"14.79,\n"
"15.2596,\n"
"15.5625,\n"
"15.9033,\n"
"16.0921,\n"
"16.3906,\n"
"16.7547,\n"
"16.9901,\n"
"17.1289,\n"
"17.3856,\n"
"17.7084,\n"
"17.9363,\n"
"17.9533,\n"
"18.1325,\n"
"18.3552,\n"
"18.4859,\n"
"18.5357,\n"
"18.5249,\n"
"18.5048,\n"
"18.622,\n"
"18.4414,\n"
"18.6183,\n"
"18.5074,\n"
"18.4496,\n"
"18.1623,\n"
"17.8825,\n"
"17.6903,\n"
"17.4322,\n"
"16.8687,\n"
"16.6486,\n"
"15.5867,\n"
"15.4747,\n"
"14.7538,\n"
"13.8882,\n"
"0.0 ] # Last entry ([101]) is 0 for kludge for high ratios.\n"
"# negativeFermionThermalFunctionValues is a list of length 102 (entries\n"
"# from [0] to [101]) with values for negative m^2/T^2 ratios for\n"
"# J_{+}( m^2/T^2 ) calculated with Mathematica.\n"
"negativeFermionThermalFunctionValues = [ fermionLowA,\n"
"-2.07364,\n"
"-2.27358,\n"
"-2.48534,\n"
"-2.70391,\n"
"-2.92541,\n"
"-3.14658,\n"
"-3.36491,\n"
"-3.57624,\n"
"-3.77925,\n"
"-3.97078,\n"
"-4.14835,\n"
"-4.30904,\n"
"-4.45004,\n"
"-4.56835,\n"
"-4.66085,\n"
"-4.7235,\n"
"-4.75264,\n"
"-4.74345,\n"
"-4.69105,\n"
"-4.58813,\n"
"-4.4273,\n"
"-4.19694,\n"
"-3.88134,\n"
"-3.45309,\n"
"-2.83773,\n"
"-1.98787,\n"
"-1.16561,\n"
"-0.36077,\n"
"0.406634,\n"
"1.15178,\n"
"1.87767,\n"
"2.571,\n"
"3.2662,\n"
"3.89582,\n"
"4.52003,\n"
"5.12808,\n"
"5.70342,\n"
"6.27452,\n"
"6.80631,\n"
"7.34613,\n"
"7.80949,\n"
"8.27843,\n"
"8.73792,\n"
"9.20238,\n"
"9.6278,\n"
"9.97677,\n"
"10.2511,\n"
"10.6802,\n"
"11.0228,\n"
"11.381,\n"
"11.6809,\n"
"11.9586,\n"
"12.185,\n"
"12.4718,\n"
"12.6165,\n"
"12.8352,\n"
"13.0107,\n"
"13.2081,\n"
"13.4074,\n"
"13.5214,\n"
"13.68,\n"
"13.7295,\n"
"13.844,\n"
"13.9108,\n"
"13.9615,\n"
"13.9883,\n"
"14.0893,\n"
"14.0869,\n"
"14.0305,\n"
"14.0015,\n"
"13.9733,\n"
"14.0218,\n"
"14.1443,\n"
"13.8104,\n"
"13.7344,\n"
"13.6113,\n"
"13.4642,\n"
"13.339,\n"
"13.0696,\n"
"13.2315,\n"
"13.2977,\n"
"12.8528,\n"
"13.2042,\n"
"12.7854,\n"
"12.1753,\n"
"11.9947,\n"
"11.7797,\n"
"11.8132,\n"
"11.4267,\n"
"10.993,\n"
"11.0223,\n"
"10.6826,\n"
"10.103,\n"
"9.8427,\n"
"9.46677,\n"
"9.40834,\n"
"9.15316,\n"
"8.51647,\n"
"8.1202,\n"
"8.22949,\n"
"0.0 ] # Last entry ([101]) is 0 for kludge for high ratios.\n"
"largestSquareRatioMagnitude = ( 2.0 * math.pi )**2\n"
"negativeThermalInterpolationStepWidth = ( 0.01\n"
"                                          * largestSquareRatioMagnitude )\n"
"\n"
"def NegativeThermalFunction( squareRatio, interpolationList ):\n"
"    scaledRatio = ( abs( squareRatio * 100.0 )\n"
"                    / largestSquareRatioMagnitude )\n"
"    if ( scaledRatio > 100.0 ):\n"
"        if False:\n"
"          print( \"Warning! Thermal function called with negative ratio\"\n"
"                 + \" of m^2/T^2 with very large magnitude (\"\n"
"                 + str( squareRatio )\n"
"                 + \"). The interpolation function only goes to\"\n"
"                 + \" -(2 pi)^2, and 0.0 will be returned, as the\"\n"
"                 + \" relative contribution of the thermal term for this\"\n"
"                 + \" degree of freedom should be negligible compared to\"\n"
"                 + \" the tree-level part for such a large value of\"\n"
"                 + \" |m^2/T^2|.\" )\n"
"        if ( scaledRatio > 101.0 ):\n"
"            scaledRatio = 100.9\n"
"# The indices looked for are capped to [100] and [101], which is fine as\n"
"# the lists are 102 elements long.\n"
"    ratioFloor = math.floor( scaledRatio )\n"
"    floorIndex = int( ratioFloor )\n"
"    ratioRemainder = ( scaledRatio - ratioFloor )\n"
"    constantValue = interpolationList[ floorIndex ]\n"
"    slopeValue = ( ( interpolationList[ floorIndex + 1 ] - constantValue )\n"
"                   / negativeThermalInterpolationStepWidth )\n"
"    return ( ( interpolationList[ floorIndex ] * ( 1.0 - ratioRemainder ) )\n"
"             + ( interpolationList[ floorIndex + 1 ] * ratioRemainder ) )\n"
"\n"
"\n"
"\n"
"# The function ThermalFunction returns -J_{-} or +J{+}, so should be added\n"
"# to the potential regardless of whether the degree of freedom is\n"
"# bosonic (+) or fermionic (-): i.e. for n bosons of the same mass-squared,\n"
"# add +n * ThermalFunction, while for n fermions of the same mass-squared,\n"
"# add +n * ThermalFunction.\n"
"def ThermalFunction( massSquared,\n"
"                     temperatureSquared,\n"
"                     signedDegreesOfFreedom ):\n"
"    if ( massSquared < 0.0 ):\n"
"        if ( signedDegreesOfFreedom < 0 ):\n"
"            return NegativeThermalFunction( ( massSquared\n"
"                                              / temperatureSquared ),\n"
"                                     negativeFermionThermalFunctionValues )\n"
"        else:\n"
"            return NegativeThermalFunction( ( massSquared\n"
"                                              / temperatureSquared ),\n"
"                                       negativeBosonThermalFunctionValues )\n"
"    else:\n"
"        massOverTemperature = math.sqrt( massSquared / temperatureSquared )\n"
"        if ( massOverTemperature > upperKinkRatio ):\n"
"            return HighRatioApproximation( massOverTemperature )\n"
"        else:\n"
"            if ( massOverTemperature < lowerKinkRatio ):\n"
"                if ( signedDegreesOfFreedom < 0 ):\n"
"                    return LowRatioFermionApproximation(\n"
"                                                      massOverTemperature )\n"
"                else:\n"
"                    return LowRatioBosonApproximation(\n"
"                                                      massOverTemperature )\n"
"            else:\n"
"                if ( signedDegreesOfFreedom < 0 ):\n"
"                    return IntermediateRatioFermionApproximation(\n"
"                                                      massOverTemperature )\n"
"                else:\n"
"                    return IntermediateRatioBosonApproximation(\n"
"                                                      massOverTemperature )\n"
"\n"
"\n"
"def ThermalCorrections( massSquaredValues,\n"
"                        overallFactor,\n"
"                        temperatureSquared ):\n"
"    thermalCorrection = 0.0\n"
"    for massSquaredValue in massSquaredValues:\n"
"        thermalCorrection += ThermalFunction( massSquaredValue.real,\n"
"                                              temperatureSquared,\n"
"                                              overallFactor )\n"
"    return ( abs( overallFactor ) * thermalCorrection )\n"
"\n"
"\n"
"def AbsThermalCorrections( massSquaredValues,\n"
"                           overallFactor,\n"
"                           temperatureSquared ):\n"
"    thermalCorrection = 0.0\n"
"    for massSquaredValue in massSquaredValues:\n"
"        thermalCorrection += ThermalFunction( abs( massSquaredValue.real ),\n"
"                                              temperatureSquared,\n"
"                                              overallFactor )\n"
"    return ( abs( overallFactor ) * thermalCorrection )\n"
"\n"
"\n"
"def FloorThermalCorrections( massSquaredValues,\n"
"                             overallFactor,\n"
"                             temperatureSquared ):\n"
"    thermalCorrection = 0.0\n"
"    for massSquaredValue in massSquaredValues:\n"
"        massSquared = massSquaredValue.real\n"
"        if ( massSquared < 0.0 ):\n"
"            massSquared = 0.0\n"
"        thermalCorrection += ThermalFunction( massSquared,\n"
"                                              temperatureSquared,\n"
"                                              overallFactor )\n"
"    return ( abs( overallFactor ) * thermalCorrection )\n"
"\n"
"\n"
"thermalFactor = ( 1.0 / ( 2.0 * math.pi * math.pi ) )\n"
"\n"
<< potentialMinimizer.prepareParameterDependentPython( sarahInterpreter )
<< "\n"
"# Vevacious class:\n"
"\n"
"class Vevacious:\n"
"    \"\"\"\n"
"    This class holds all the parameter-dependent stuff such as various\n"
"    refinements of the potential energy as a function of field\n"
"    configurations, as well as ancillary functions.\n"
"    \"\"\"\n"
"\n"
"    def __init__( self,\n"
"                  EffectivePotential,\n"
"                  ageOfKnownUniverseInSeconds = 4.3E+17,\n"
"                  fourthRootOfSolitonicFactorA = energyScale,\n"
"                  tunnelPathResolutionInGev = 50.0,\n"
"                  currentTemperature = 0.0 ):\n"
"        self.vevaciousVersion = \"" << vevaciousVersion << "\"\n"
"        self.outputFile = \"./VevaciousResult.vout\"\n"
"        self.warningMessages = []\n"
"\n"
"        self.EffectivePotential = EffectivePotential\n"
"        self.potentialScaler = PotentialScaler( self.EffectivePotential,\n"
"                                    temperatureValue = currentTemperature )\n"
"        self.PyMinuit = minuit.Minuit(\n"
"                   self.potentialScaler.ScaledFunctionFromScaledArguments )\n"
"        self.SetFieldValueLimit( 10.0 * energyScale )\n"
"        self.rollingToleranceSquared = ( " << rollingTolerance << " )**2\n"
"        self.foundMinima = []\n"
"        self.saddleSplitNudges = [ ";
    for( std::vector< double >::iterator
         whichNudge( saddleNudgeList.begin() );
         saddleNudgeList.end() > whichNudge;
         ++whichNudge )
    {
      if( saddleNudgeList.begin() != whichNudge )
      {
        outputFile << ", ";
      }
      outputFile << "( " << *whichNudge << " )";
    }
    outputFile << " ]\n"
"        self.fieldScaling = 2.0\n"
"        self.dsbVacuumIsMetastable = False\n"
"        self.originAsExtremum = self.FieldDictionaryToExtremum(\n"
"                                                              fieldOrigin )\n"
"        self.dsbVacuum = self.TryToMinimize( dsbInput )\n"
"        self.dsbArray = FieldDictionaryToArray( self.dsbVacuum[\n"
"                                                        \"FieldValues\" ] )\n"
"        self.minimumSeparationSquared = ( self.rollingToleranceSquared\n"
"                                        * numpy.sum( (self.dsbArray)**2 ) )\n"
"        self.panicVacuum = self.dsbVacuum\n"
"        self.panicArray = self.dsbArray\n"
"        self.globalMinimum = self.dsbVacuum\n"
"        self.globalArray = self.dsbArray\n"
"\n"
"        self.ageOfKnownUniverseInSeconds = ageOfKnownUniverseInSeconds\n"
"        self.ageOfKnownUniverseInInverseGev = ( ageOfKnownUniverseInSeconds\n"
"                                                / 6.58211928E-25 )\n"
"        self.fourthRootOfSolitonicFactorA = fourthRootOfSolitonicFactorA\n"
"        self.tunnelPathResolutionInGev = tunnelPathResolutionInGev\n"
"        self.currentQuantumAction = None\n"
"        self.previousQuantumAction = None\n"
"        self.currentThermalAction = None\n"
"        self.previousThermalAction = None\n"
"        self.tunnelingTimeBound = " << lifetimeThreshold << "\n"
"        self.quantumActionThreshold = -1.0\n"
"        if ( 0.0 < self.tunnelingTimeBound ):\n"
"            self.quantumActionThreshold = ( 4.0\n"
"                                        * math.log( self.tunnelingTimeBound\n"
"                                      * self.ageOfKnownUniverseInInverseGev\n"
"                                         * fourthRootOfSolitonicFactorA ) )\n"
"        self.thermalProbabilityThreshold = " << survivalThreshold << "\n"
"# Based on correspondence with Alexander Kusenko and discussion with Bjoern\n"
"# Garbrecht:\n"
"# Taking [decay width per horizon]\n"
"# = [horizon volume] * [solitonic coefficient] * exp(-[thermal action]/T)\n"
"# at temperature T, where [horizon volume] = ( M_Plank / T^2 )^3, and\n"
"# taking [solitonic coefficient] to be T^4,\n"
"# the survival probability per horizon =\n"
"# exp( -integral of [time at T] with respect to [decay time] )\n"
"# = exp( -integral of [decay width per horizon] dT * [factor] ) )\n"
"# which exponents for N horizons to\n"
"# exp( -N * integral * [factor] ) )\n"
"# and [decay width per horizon] = M_Plank^3 T^(-2) exp(-S_3(T)/T)\n"
"# where [thermal action at temperature T] = S_3(T)\n"
"# exp( -N * integral * [factor] ) ) can be written, from entropy\n"
"# conservation and so on, as\n"
"# exp( -N * integral of C T^(-2) exp(-S_3(T)/T) dT ) )\n"
"# where C = [reduced Planck mass] * [solitonic coefficient/T^4]\n"
"# * sqrt[45/(4 pi^3 g_star(T))] * [g_star^now/g_star(T)] * (T_now/H_now)^3\n"
"# and we take g_star(T) to be 105.75 (what it is for the SM above\n"
"# temperatures of m_top) and conservatively take it as constant from T = 0\n"
"# to T_opt. Hence we have\n"
"# exp( -1.581E+106 GeV * integral of T^(-2) exp(-S_3(T)/T) dT ) )\n"
"# integrated from T = 0 to T_opt (as the contribution from higher\n"
"# temperatures drops off very quickly).\n"
"# 1.581E+106 is exp( 244.53 = lnOfThermalIntegrationFactor ) which is\n"
"# in agreement with the value of 240 quoted in the CosmoTransitions manual\n"
"# for an estimate of the threshold S_3(T)/T for T= 100 GeV.\n"
"# Kusenko (and others in the literature, including Wainwright implicitly in\n"
"# the CosmoTransitions manual as just mentioned) took the integral of\n"
"# exp(-S_3(T)/T) T^(-2) to be exp( S_3(T_opt)/T_opt) T_opt^(-1) where T_opt\n"
"# is the optimal tunneling temperature which dominates the integral.\n"
"# This might be a bit aggressive, and taking S_3(T) to be approximated by\n"
"# S_3(0) + T S' leads to the integral being \n"
"# exp( -S_3(T_opt)/T_opt ) / S_3(0)\n"
"# Assuming that S_3(0) < S_3(T_opt) (which should hold for all cases of\n"
"# interest), the full integral should be between\n"
"# exp( -S_3(T_opt)/T_opt ) / S_3(T_opt)\n"
"# and exp( -S_3(T_opt)/T_opt ) / T_opt\n"
"# For a threshold survival probability P,\n"
"# 1.581E+106 GeV * integral of T^(-2) exp(-S_3(T)/T) dT )\n"
"# should be larger than ln(1/P). Hence we compare\n"
"# (S_3(T_opt)/T_opt) + ln( x / GeV ) to\n"
"# lnOfThermalIntegrationFactor - ln( ln(1/P) ) where x is either\n"
"# S_3(T_opt) or T_opt.\n"
"        self.thermalActionOverTemperatureComparison = (\n"
"                                               lnOfThermalIntegrationFactor\n"
"                                                                - math.log(\n"
"                                                              math.log( 1.0\n"
"                                   / self.thermalProbabilityThreshold ) ) )\n"
"        self.startingTime = time.clock()\n"
"        self.allowedRunningTime = 100.0\n"
"        self.ShouldTunnel = " << shouldTunnel << "\n"
"        self.ShouldTunnelThermally = " << tunnelThermally << "\n"
"        self.ShouldDeformTunnelPaths = " << deformTunnelPaths << "\n"
"        if ( self.minimumSeparationSquared < 1.0 ):\n"
"            self.minimumSeparationSquared = 1.0\n"
"\n"
"\n"
"    def SetFieldValueLimit( self, fieldValueLimit ):\n"
"        for fieldVariable in self.PyMinuit.values:\n"
"            self.PyMinuit.limits[ fieldVariable ] = ( -fieldValueLimit,\n"
"                                                      fieldValueLimit )\n"
"\n"
"\n"
"    def MinuitValuesFromDictionary( self, fieldValueDictionary ):\n"
"        scaledFields = fieldValueDictionary.copy()\n"
"        for fieldKey in scaledFields.keys():\n"
"            scaledFields[ fieldKey ] *= inverseScale\n"
"        return scaledFields\n"
"\n"
"\n"
"    def FieldValuesFromMinuit( self, fieldValueDictionary ):\n"
"        scaledFields = fieldValueDictionary.copy()\n"
"        for fieldKey in scaledFields.keys():\n"
"            scaledFields[ fieldKey ] *= energyScale\n"
"        return scaledFields\n"
"\n"
"\n"
"    def PotentialValueFromMinuit( self, potentialDepth ):\n"
"        return ( potentialDepth * energyScaleFourth )\n"
"\n"
"\n"
"    def FieldDictionaryToExtremum( self, fieldValueDictionary ):\n"
"        return { \"FieldValues\": fieldValueDictionary.copy(),\n"
"                 \"PotentialValue\": FunctionFromDictionary(\n"
"                                                   self.EffectivePotential,\n"
"                                                    fieldValueDictionary ),\n"
"                 \"MinuitError\": None }\n"
"\n"
"\n"
"    def ExtremumAsMathematica( self, givenExtremum ):\n"
"        return ( \"{ \"\n"
"             + UserFieldsAsMathematica( givenExtremum[ \"FieldValues\" ]  )\n"
"                 + \", TreeLevelPotentialValue -> \"\n"
"                 + str( FunctionFromDictionary( TreeLevelPotential,\n"
"                                       givenExtremum[ \"FieldValues\" ] ) )\n"
"                 + \", EffectivePotentialValue -> \"\n"
"                 + str( givenExtremum[ \"PotentialValue\" ] )\n"
"                 + \" }\" )\n"
"\n"
"\n"
"    def TryToMinimize( self, fieldValueDictionary ):\n"
"        foundExtremum = self.FieldDictionaryToExtremum(\n"
"                                                     fieldValueDictionary )\n"
"        print( \"trying to minimize \"\n"
"               + self.ExtremumAsMathematica( foundExtremum )\n"
"               + \" (T = \" + str( self.potentialScaler.temperatureValue )\n"
"               + \")\" )\n"
"        try:\n"
"            self.PyMinuit.values = self.MinuitValuesFromDictionary(\n"
"                                                     fieldValueDictionary )\n"
"            self.PyMinuit.migrad()\n"
"            foundExtremum = { \"FieldValues\": self.FieldValuesFromMinuit(\n"
"                                                    self.PyMinuit.values ),\n"
"                         \"PotentialValue\": self.PotentialValueFromMinuit(\n"
"                                                      self.PyMinuit.fval ),\n"
"                            \"MinuitError\": self.PotentialValueFromMinuit(\n"
"                                                      self.PyMinuit.edm ) }\n"
"            print( \"rolled to \"\n"
"                   + self.ExtremumAsMathematica( foundExtremum )\n"
"                   + \"\\n\\n\" )\n"
"        except minuit.MinuitError as minuitError:\n"
"            self.LogWarning( \"PyMinuit had problems starting at \"\n"
"                        + UserFieldsAsMathematica( fieldValueDictionary )\n"
"                        + \"! [minuit.MinuitError: \"\n"
"                        + str( minuitError )\n"
"                        + \"]. PyMinuit stopped at \"\n"
"                        + UserFieldsAsMathematica(\n"
"                       self.FieldValuesFromMinuit( self.PyMinuit.values ) )\n"
"                        + \".  Minuit\'s estimate of how much deeper it\"\n"
"                        + \" should go is \"\n"
"                + str( self.PotentialValueFromMinuit( self.PyMinuit.edm ) )\n"
"                        + \" GeV^4.\" )\n"
"        except Exception as unexpectedException:\n"
"            self.LogWarning( \"Some lower-level exception happened:\\n\"\n"
"                        + str( unexpectedException )\n"
"                        + \"\\nwhile PyMinuit was trying to minimize the\"\n"
"                        + \" potential starting at \"\n"
"                        + UserFieldsAsMathematica( fieldValueDictionary )\n"
"                        + \". PyMinuit stopped at \"\n"
"                        + UserFieldsAsMathematica(\n"
"                       self.FieldValuesFromMinuit( self.PyMinuit.values ) )\n"
"                        + \".  Minuit\'s estimate of how much deeper it\"\n"
"                        + \" should go is \"\n"
"                + str( self.PotentialValueFromMinuit( self.PyMinuit.edm ) )\n"
"                        + \" GeV^4.\" )\n"
"        return foundExtremum\n"
"\n"
"\n"
"    def LogWarning( self, warningMessage ):\n"
"        self.warningMessages.append( warningMessage )\n"
"        print( warningMessage )\n"
"\n"
"\n"
"    def WriteExtrema( self, fieldValueDictionaries, fileToWrite ):\n"
"        outputFile = open( fileToWrite, \"w\" )\n"
"        outputFile.write( \"{\\n\" )\n"
"        for fieldValueDictionary in fieldValueDictionaries:\n"
"            currentExtremum = self.FieldDictionaryToExtremum(\n"
"                                                     fieldValueDictionary )\n"
"            outputFile.write( self.ExtremumAsMathematica( currentExtremum )\n"
"                              + \"\\n\" )\n"
"        outputFile.write( \"}\" )\n"
"        outputFile.close()\n"
"\n"
"\n"
"    def CloseEnoughToDsbOrSignFlip( self, fieldValueDictionary ):\n"
"        fieldArray = FieldDictionaryToArray( fieldValueDictionary )\n"
"        longerLengthSquared = numpy.sum( fieldArray**2 )\n"
"        otherLengthSquared = numpy.sum( (self.dsbArray)**2 )\n"
"        if( otherLengthSquared > longerLengthSquared ):\n"
"            longerLengthSquared = otherLengthSquared\n"
"        thresholdLengthSquared = ( longerLengthSquared\n"
"                                   * self.rollingToleranceSquared )\n"
"        return ( numpy.sum( ( abs( fieldArray )\n"
"                              - abs( self.dsbArray ) )**2 )\n"
"                 < thresholdLengthSquared )\n"
"\n"
"\n"
"    def TryToMinimizeIncludingRescaling( self, fieldValueDictionary ):\n"
"        foundExtremum = self.TryToMinimize( fieldValueDictionary )\n"
"        if ( self.CloseEnoughToDsbOrSignFlip( foundExtremum[\n"
"                                                        \"FieldValues\" ] )\n"
"             and\n"
"             not self.CloseEnoughToDsbOrSignFlip( fieldValueDictionary ) ):\n"
"            scaledValues = fieldValueDictionary.copy()\n"
"            for fieldKey in fieldValueDictionary.keys():\n"
"                scaledValues[ fieldKey ] *= self.fieldScaling\n"
"            scaledPoint = self.FieldDictionaryToExtremum( scaledValues )\n"
"            self.LogWarning( \"Starting point differed from DSB vacuum\"\n"
"                             + \" but rolled there. Trying again with\"\n"
"                             + \" scaled field configuration (\"\n"
"                             + self.ExtremumAsMathematica(\n"
"                                            self.FieldDictionaryToExtremum(\n"
"                                         scaledPoint[ \"FieldValues\" ] ) )\n"
"                             + \").\" )\n"
"            foundExtremum = self.TryToMinimize( scaledValues )\n"
"        return foundExtremum\n"
"\n"
"\n"
"    def NudgeOff( self, foundSaddlePoint, saddleSplitNudge ):\n"
"        saddleFieldValues = foundSaddlePoint[ 0 ][ \"FieldValues\" ]\n"
"        positivelyNudged = saddleFieldValues.copy()\n"
"        negativelyNudged = saddleFieldValues.copy()\n"
"        nudgeVector = foundSaddlePoint[ 1 ]\n"
"        for fieldIndex in range( len( namesOfFields ) ):\n"
"            whichField = namesOfFields[ fieldIndex ]\n"
"            positivelyNudged[ whichField ] += nudgeVector[ fieldIndex ]\n"
"            negativelyNudged[ whichField ] -= nudgeVector[ fieldIndex ]\n"
"        return [ positivelyNudged, negativelyNudged ]\n"
"\n"
"\n"
"    def RollExtrema( self, startingFieldConfigurations ):\n"
"        foundExtrema = []\n"
"        foundSaddlePoints = []\n"
"        for startingFieldConfiguration in startingFieldConfigurations:\n"
"            foundExtrema.append( self.TryToMinimizeIncludingRescaling(\n"
"                                             startingFieldConfiguration ) )\n"
"        if ( len( self.saddleSplitNudges ) > 0 ):\n"
"            for foundExtremum in foundExtrema:\n"
"                isMinimum, nudgeVector = self.CheckHessian( foundExtremum )\n"
"                if ( isMinimum ):\n"
"                    self.foundMinima.append( foundExtremum )\n"
"                else:\n"
"                    foundSaddlePoints.append( [ foundExtremum,\n"
"                                                nudgeVector ] )\n"
"        else:\n"
"            self.foundMinima = foundExtrema\n"
"# At this point, the field configurations from startingFieldConfigurations\n"
"# have been divided into minima in self.foundMinima and saddle points (&\n"
"# possibly maxima) in foundSaddlePoints, if self.saddleSplitNudges is not\n"
"# an empty list. (If self.saddleSplitNudges is an empty list, at this point\n"
"# self.foundMinima contains all the extrema from the 1st pass, & this\n"
"# function won\'t do anything else.)\n"
"        for saddleSplitNudge in self.saddleSplitNudges:\n"
"# For each requested nudging, nudgedFieldConfigurations is made from all\n"
"# the pairs of nearby configurations to each in foundSaddlePoints,\n"
"# foundSaddlePoints is then emptied, & each configuration is rolled & the\n"
"# rolled configuration is put into either self.foundMinima if a minimum or\n"
"# foundSaddlePoints again if not, & then the loop iterates.\n"
"            nudgedFieldConfigurations = []\n"
"            for foundSaddlePoint in foundSaddlePoints:\n"
"                nudgedFieldConfigurations.extend( self.NudgeOff(\n"
"                                                          foundSaddlePoint,\n"
"                                                       saddleSplitNudge ) )\n"
"            foundSaddlePoints = []\n"
"            for nudgedFieldConfiguration in nudgedFieldConfigurations:\n"
"                foundExtremum = self.TryToMinimizeIncludingRescaling(\n"
"                                                 nudgedFieldConfiguration )\n"
"                isMinimum, nudgeVector = self.CheckHessian( foundExtremum )\n"
"                if ( isMinimum ):\n"
"                    self.foundMinima.append( foundExtremum )\n"
"                else:\n"
"                    foundSaddlePoints.append( [ foundExtremum,\n"
"                                                nudgeVector ] )\n"
"# At this point, any remaining saddle points / maxima are just appended to\n"
"# self.foundMinima so that they can still be tunneled to, regardless of the\n"
"# fact that the Universe would still evolve away from such configurations.\n"
"        for foundSaddlePoint in foundSaddlePoints:\n"
"            self.foundMinima.append( foundSaddlePoint[ 0 ] )\n"
"# Now the minima are sorted by depth:\n"
"        tunnelDistance = -1.0\n"
"        for foundMinimum in self.foundMinima:\n"
"            if ( ( foundMinimum[ \"PotentialValue\" ] < self.dsbVacuum[\n"
"                                                    \"PotentialValue\" ]  )\n"
"                 and\n"
"                 ( not self.CloseEnoughToDsbOrSignFlip( foundMinimum[\n"
"                                                   \"FieldValues\" ] ) ) ):\n"
"                self.dsbVacuumIsMetastable = True\n"
"                fieldArray = FieldDictionaryToArray( foundMinimum[\n"
"                                                        \"FieldValues\" ] )\n"
"                currentDistanceToDsb = numpy.sum( ( fieldArray\n"
"                                                    - self.dsbArray )**2 )\n"
"                if ( ( tunnelDistance < 0.0 )\n"
"                     or\n"
"                     ( currentDistanceToDsb < tunnelDistance ) ):\n"
"                    tunnelDistance = currentDistanceToDsb\n"
"                    self.panicVacuum = foundMinimum\n"
"                    self.panicArray = fieldArray\n"
"                if ( foundMinimum[ \"PotentialValue\"\n"
"                            ] < self.globalMinimum[ \"PotentialValue\" ] ):\n"
"                    self.globalMinimum = foundMinimum\n"
"                    self.globalArray = fieldArray\n"
"\n"
"\n"
"    def CheckHessian( self, foundExtremum ):\n"
"        isMinimum = True\n"
"        nudgeVector = list( originArray )\n"
"        eigensystemOfHessian = numpy.linalg.eigh( NumericalHessian(\n"
"                                                   self.EffectivePotential,\n"
"                                       foundExtremum[ \"FieldValues\" ] ) )\n"
"# It doesn\'t matter what value mostNegativeEigenvalueValue has if it\'s\n"
"# positive because it is only there to record the lowest eigenvalue if it\n"
"# is less than or equal to zero.\n"
"        mostNegativeEigenvalueValue = 1.0\n"
"        mostNegativeEigenvalueIndex = 0\n"
"        for eigenvalueIndex in range( len( eigensystemOfHessian[ 0 ] ) ):\n"
"            if ( mostNegativeEigenvalueValue > eigensystemOfHessian[ 0 ][\n"
"                                                       eigenvalueIndex ] ):\n"
"                mostNegativeEigenvalueValue = eigensystemOfHessian[ 0 ][\n"
"                                                          eigenvalueIndex ]\n"
"                mostNegativeEigenvalueIndex = eigenvalueIndex\n"
"        isMinimum = ( mostNegativeEigenvalueValue > 0.0 )\n"
"        if not isMinimum:\n"
"            nudgeVector = eigensystemOfHessian[ 1 ][ :,\n"
"                                              mostNegativeEigenvalueIndex ]\n"
"        return [ isMinimum, nudgeVector ]\n"
"\n"
"\n"
"    def WriteMinima( self, fileToWrite ):\n"
"        outputFile = open( fileToWrite, \"w\" )\n"
"        outputFile.write( \"{\\n\" )\n"
"        for foundMinimum in self.foundMinima:\n"
"            outputFile.write( self.ExtremumAsMathematica( foundMinimum )\n"
"                              + \"\\n\" )\n"
"        outputFile.write( \"}\" )\n"
"        outputFile.close()\n"
"\n"
"\n"
"    def CalculateAction( self,\n"
"                         falseVacuum = None,\n"
"                         trueVacuum = None,\n"
"                         deformPath = True,\n"
"                         thermalNotQuantum = None,\n"
"                         innerLoopMaxDeformations = 10,\n"
"                         outerLoopMaxDeformations = 10 ):\n"
"        if ( falseVacuum is None ):\n"
"            falseVacuum = self.dsbVacuum\n"
"        if ( trueVacuum is None ):\n"
"            trueVacuum = self.panicVacuum\n"
"        if thermalNotQuantum is None:\n"
"            thermalNotQuantum = (\n"
"                            self.potentialScaler.currentTemperature > 0.0 )\n"
"        falseVacuumArray = FieldDictionaryToArray( falseVacuum[\n"
"                                                        \"FieldValues\" ] )\n"
"        trueVacuumArray = FieldDictionaryToArray( trueVacuum[\n"
"                                                        \"FieldValues\" ] )\n"
"        arrayOfArrays = numpy.array( [ trueVacuumArray,\n"
"                                                       falseVacuumArray ] )\n"
"        configurationDifferenceLength = math.sqrt( numpy.sum( (\n"
"                                                            trueVacuumArray\n"
"                                                - falseVacuumArray )**2 ) )\n"
"        tunnelPathPoints = int( configurationDifferenceLength\n"
"                                / self.tunnelPathResolutionInGev )\n"
"        print( \"tunnelingCalculator.npoints being given \"\n"
"               + str( tunnelPathPoints ) )\n"
"        print( \"falseVacuumArray = \" + str( falseVacuumArray ) )\n"
"        print( \"trueVacuumArray = \" + str( trueVacuumArray ) )\n"
"\n"
"        numberOfSpaceTimeDimensionsForTunneling = 3\n"
"\n"
"        if thermalNotQuantum:\n"
"            self.previousThermalAction = self.currentThermalAction\n"
"            numberOfSpaceTimeDimensionsForTunneling = 2\n"
"        else:\n"
"            self.previousQuantumAction = self.currentQuantumAction\n"
"        calculatedAction = exponentCutOff\n"
"\n"
"        if ( ctMajorVersion is 2 ):\n"
"            if not deformPath:\n"
"                innerLoopMaxDeformations = 0\n"
"                outerLoopMaxDeformations = 1\n"
"            tunnelingResult = CTPD.fullTunneling( path_pts = arrayOfArrays,\n"
"                              V = self.potentialScaler.PotentialFromMatrix,\n"
"                              dV = self.potentialScaler.GradientFromMatrix,\n"
"                                              tunneling_init_params = dict(\n"
"                         alpha = numberOfSpaceTimeDimensionsForTunneling ),\n"
"                                       tunneling_findProfile_params = dict(\n"
"                                              npoints = tunnelPathPoints ),\n"
"                                          deformation_deform_params = dict(\n"
"                                      maxiter = innerLoopMaxDeformations ),\n"
"                                       maxiter = outerLoopMaxDeformations )\n"
"            calculatedAction = tunnelingResult.action\n"
"        elif ( ctMajorVersion is 1 ):\n"
"            tunnelingCalculator = CTPD.fullTunneling( phi = arrayOfArrays,\n"
"                              V = self.potentialScaler.PotentialFromMatrix,\n"
"                              dV = self.potentialScaler.GradientFromMatrix,\n"
"                           alpha = numberOfSpaceTimeDimensionsForTunneling,\n"
"                                                npoints = tunnelPathPoints,\n"
"                                                   quickTunneling = False )\n"
"            if deformPath:\n"
"                tunnelingCalculator.run(\n"
"                                        maxiter = innerLoopMaxDeformations,\n"
"                                      maxiter2 = outerLoopMaxDeformations )\n"
"            else:\n"
"                tunnelingCalculator.tunnel1D()\n"
"            calculatedAction = tunnelingCalculator.findAction()\n"
"\n"
"        if thermalNotQuantum:\n"
"            self.currentThermalAction = calculatedAction\n"
"        else:\n"
"            self.currentQuantumAction = calculatedAction\n"
"\n"
"        return calculatedAction\n"
"\n"
"\n"
"    def AllowedRunningTimeExceeded( self ):\n"
"        returnValue = ( time.clock() > ( self.startingTime\n"
"                                         + self.allowedRunningTime ) )\n"
"        print( \"AllowedRunningTimeExceeded about to return \"\n"
"               + str( returnValue )\n"
"               + \".\" )\n"
"        return returnValue\n"
"\n"
"\n"
"    def QuantumTunnelingTimeInInverseGev( self, currentAction ):\n"
"        if ( currentAction < 0.0 ):\n"
"            print( \"Warning! Quantum tunneling time for a negative\"\n"
"                   + \" action (\"\n"
"                   + str( currentAction )\n"
"                   + \") asked for, which may indicate a problem (or may\"\n"
"                   + \" not, but it would be a very thin potential\"\n"
"                   + \" barrier). Returning a tunneling time for an\"\n"
"                   + \" action of 0.\" )\n"
"            currentAction = 0.0\n"
"        elif ( currentAction > exponentCutOff ):\n"
"            print( \"Warning! Quantum tunneling time for a very large\"\n"
"                   + \" action (\"\n"
"                   + str( currentAction )\n"
"                   + \") asked for, which would probably cause an\"\n"
"                   + \" overflow in the evaluation of the exponent.\"\n"
"                   + \" Returning a tunneling time for an action of \"\n"
"                   + str( exponentCutOff )\n"
"                   + \".\" )\n"
"            currentAction = exponentCutOff\n"
"        return ( math.exp( 0.25 * currentAction )\n"
"                 / self.fourthRootOfSolitonicFactorA )\n"
"\n"
"\n"
"    def SetTemperature( self, currentTemperature ):\n"
"        self.potentialScaler.temperatureValue = currentTemperature\n"
"\n"
"\n"
"    def RollVacuaAtTemperature( self,\n"
"                                temperatureGuess,\n"
"                                dsbEvaporationTemperature ):\n"
"        print( \"Trying \" + str( temperatureGuess ) + \" GeV.\\n\" )\n"
"        self.SetTemperature( temperatureGuess )\n"
"        thermalDsbVacuum = self.originAsExtremum\n"
"        if ( temperatureGuess > dsbEvaporationTemperature ):\n"
"            print( \"DSB vacuum evaporated!\" )\n"
"        else:\n"
"            print( \"DSB vacuum:\" )\n"
"            thermalDsbVacuum = self.TryToMinimize( self.dsbVacuum[\n"
"                                                        \"FieldValues\" ] )\n"
"        print( \"Panic vacuum:\" )\n"
"        thermalPanicVacuum = self.TryToMinimize( self.panicVacuum[\n"
"                                                        \"FieldValues\" ] )\n"
"        return [ thermalDsbVacuum, thermalPanicVacuum ]\n"
"\n"
"\n"
"    def TunnelingPossible( self, falseVacuum, trueVacuum ):\n"
"        if ( falseVacuum[ \"PotentialValue\" ]\n"
"             > trueVacuum[ \"PotentialValue\" ] ):\n"
"            trueArray = FieldDictionaryToArray( trueVacuum[\n"
"                                                        \"FieldValues\" ] )\n"
"            if ( numpy.sum( trueArray**2 )\n"
"                 > self.minimumSeparationSquared ):\n"
"                falseArray = FieldDictionaryToArray( falseVacuum[\n"
"                                                        \"FieldValues\" ] )\n"
"                return ( numpy.sum( ( trueArray - falseArray )**2 )\n"
"                         > self.minimumSeparationSquared )\n"
"        return False\n"
"\n"/*
"\n"
"# This estimates at what temperature tunneling from the DSB vacuum to the\n"
"# panic vacuum becomes impossible.\n"
"    def FindCriticalTemperature( self, logTwoAccuracy = 7 ):\n"
"# As a 1st guess, we take the temperature such that the SM thermal\n"
"# corrections at zero VEVs would equal the difference:\n"
"        potentialDifference = ( self.dsbVacuum[  \"PotentialValue\" ]\n"
"                                - self.panicVacuum[  \"PotentialValue\" ] )\n"
"# The corrections are ( T^4 / ( 2 pi^2 ) ) * sum of J functions, and the\n"
"# values of the J functions are about 2 for massless bosonic & fermionic\n"
"# degrees of freedom, & there are ~100 degrees of freedom in the SM. Hence\n"
"# we take the coefficient of T^4 to be 100 / ( 2 pi^2 ) ~= 5\n"
"        criticalTemperatureGuess = math.pow( ( potentialDifference / 5.0 ),\n"
"                                             0.25 )\n"
"# We aim to have a pair of temperatures, one above the critical\n"
"# temperature, the other below. If the initial guess was below the critical\n"
"# temperature, we start doubling the temperature, recording the previous\n"
"# temperature each time. If it was above, we start halving the temperature,\n"
"# recording the previous temperature each time.\n"
"        print( \"Searching for critical temperature.\\n\\n\" )\n"
"        thermalDsbVacuum, thermalPanicVacuum = self.RollVacuaAtTemperature(\n"
"                                                 criticalTemperatureGuess )\n"
"        while ( not self.TunnelingPossible( thermalDsbVacuum,\n"
"                                            thermalPanicVacuum ) ):\n"
"            criticalTemperatureGuess = ( 0.5 * criticalTemperatureGuess )\n"
"            [ thermalDsbVacuum,\n"
"              thermalPanicVacuum ] = self.RollVacuaAtTemperature(\n"
"                                                 criticalTemperatureGuess )\n"
"# At this point, criticalTemperatureGuess is definitely below the critical\n"
"# temperature. Now we look for the critical temperature from below:\n"
"        while self.TunnelingPossible( thermalDsbVacuum,\n"
"                                      thermalPanicVacuum ):\n"
"            criticalTemperatureGuess = ( 2.0 * criticalTemperatureGuess )\n"
"            [ thermalDsbVacuum,\n"
"              thermalPanicVacuum ] = self.RollVacuaAtTemperature(\n"
"                                                 criticalTemperatureGuess )\n"
"# At this point, criticalTemperatureGuess should be between 1.0 and 2.0\n"
"# times the critical temperature.\n"
"        justBelowCriticalTemperature = ( 0.5 * criticalTemperatureGuess )\n"
"        justAboveCriticalTemperature = criticalTemperatureGuess\n"
"        for narrowingStep in range( logTwoAccuracy ):\n"
"            criticalTemperatureGuess = math.sqrt(\n"
"                                               justBelowCriticalTemperature\n"
"                                           * justAboveCriticalTemperature )\n"
"            [ thermalDsbVacuum,\n"
"              thermalPanicVacuum ] = self.RollVacuaAtTemperature(\n"
"                                                 criticalTemperatureGuess )\n"
"            if self.TunnelingPossible( thermalDsbVacuum,\n"
"                                       thermalPanicVacuum ):\n"
"                justBelowCriticalTemperature = criticalTemperatureGuess\n"
"            else:\n"
"                justAboveCriticalTemperature = criticalTemperatureGuess\n"
"        print( \"Critical temperature lies between \"\n"
"               + str( justBelowCriticalTemperature )\n"
"               + \" GeV and \"\n"
"               + str( justAboveCriticalTemperature )\n"
"               + \" GeV.\" )\n"
"        return [ justBelowCriticalTemperature,\n"
"                 justAboveCriticalTemperature ]\n"
"\n"*/
"\n"
"# This estimates at what temperature tunneling from the field origin to the\n"
"# given vacuum becomes impossible. It returns a temperature just below,\n"
"# a temperature just above. A pair of zeroes are returned if the given\n"
"# minimum is not deeper than the field origin at zero temperature. It is\n"
"# not technically the evaporation temperature necessarily, as the minimum\n"
"# may still exist at this temperature, though it is not deeper than the\n"
"# purely symmetric phase vacuum, so it is unreasonable to think that the\n"
"# Universe would cool into the false vacuum at that temperature rather\n"
"# than into the symmetric minimum.\n"
"    def FindEvaporationTemperature( self,\n"
"                                    minimumAtZeroTemperature,\n"
"                                    logTwoAccuracy = 7,\n"
"                   criticalComparedToOriginRatherThanEvaporation = False ):\n"
"# As a 1st guess, we take the temperature such that the SM thermal\n"
"# corrections at zero VEVs would equal the difference (remembering that the\n"
"# value of the potential at the origin has already been subtracted from\n"
"# minimumAtZeroTemperature[ \"PotentialValue\" ]):\n"
"        zeroTemperatureDepth = minimumAtZeroTemperature[\n"
"                                                       \"PotentialValue\" ]\n"
"        if ( zeroTemperatureDepth >= 0.0 ):\n"
"            return [ 0.0, 0.0 ]\n"
"# The corrections are ( T^4 / ( 2 pi^2 ) ) * sum of J functions, and the\n"
"# values of the J functions are about 2 for massless bosonic & fermionic\n"
"# degrees of freedom, & there are ~100 degrees of freedom in the SM. Hence\n"
"# we take the coefficient of T^4 to be 100 / ( 2 pi^2 ) ~= 5\n"
"        evaporationTemperatureGuess = math.pow( ( -zeroTemperatureDepth\n"
"                                                  / 5.0 ),\n"
"                                                 0.25 )\n"
"# We aim to have a pair of temperatures, one above the critical\n"
"# temperature, the other below. If the initial guess was below the critical\n"
"# temperature, we start doubling the temperature, recording the previous\n"
"# temperature each time. If it was above, we start halving the temperature,\n"
"# recording the previous temperature each time.\n"
"        originDepth = ( FunctionFromDictionary( self.EffectivePotential,\n"
"                                                fieldOrigin,\n"
"                                              evaporationTemperatureGuess )\n"
"                        - self.potentialScaler.functionAtOrigin )\n"
"        print( \"Searching for evaporation temperature.\\n\\n\" )\n"
"        print( \"Trying \"\n"
"               + str( evaporationTemperatureGuess )\n"
"               + \" GeV.\\n\" )\n"
"        while ( originDepth > zeroTemperatureDepth ):\n"
"            evaporationTemperatureGuess = ( 2.0\n"
"                                            * evaporationTemperatureGuess )\n"
"            print( \"Trying \"\n"
"                   + str( evaporationTemperatureGuess )\n"
"                   + \" GeV.\\n\" )\n"
"            originDepth = ( FunctionFromDictionary(\n"
"                                                   self.EffectivePotential,\n"
"                                                    fieldOrigin,\n"
"                                              evaporationTemperatureGuess )\n"
"                            - self.potentialScaler.functionAtOrigin )\n"
"# At this point, criticalTemperatureGuess is definitely above the critical\n"
"# temperature. Now we look for the critical temperature from above:\n"
"        while ( originDepth < zeroTemperatureDepth ):\n"
"            evaporationTemperatureGuess = ( 0.5\n"
"                                            * evaporationTemperatureGuess )\n"
"            print( \"Trying \"\n"
"                   + str( evaporationTemperatureGuess )\n"
"                   + \" GeV.\\n\" )\n"
"            originDepth = ( FunctionFromDictionary(\n"
"                                                   self.EffectivePotential,\n"
"                                                    fieldOrigin,\n"
"                                              evaporationTemperatureGuess )\n"
"                            - self.potentialScaler.functionAtOrigin )\n"
"# At this point, evaporationTemperatureGuess should be between 1.0 and 2.0\n"
"# times the temperature at which the field origin is lower than that of the\n"
"# minimum at zero temperature.\n"
"        print( \"The field origin becomes deeper than the minimum at\"\n"
"               + \" T = 0 is between \"\n"
"               + str( evaporationTemperatureGuess )\n"
"               + \" GeV and \"\n"
"               + str( 2.0 * evaporationTemperatureGuess )\n"
"               + \" GeV.\\n\" )\n"
"        self.SetTemperature( evaporationTemperatureGuess )\n"
"        thermalMinimum = self.TryToMinimize( minimumAtZeroTemperature[\n"
"                                                        \"FieldValues\" ] )\n"
"        keepLooping = True\n"
"        while keepLooping:\n"
"            evaporationTemperatureGuess = ( 2.0\n"
"                                            * evaporationTemperatureGuess )\n"
"            print( \"Trying \"\n"
"                   + str( evaporationTemperatureGuess )\n"
"                   + \" GeV.\\n\" )\n"
"            self.SetTemperature( evaporationTemperatureGuess )\n"
"            thermalMinimum = self.TryToMinimize( minimumAtZeroTemperature[\n"
"                                                        \"FieldValues\" ] )\n"
"            if criticalComparedToOriginRatherThanEvaporation:\n"
"                keepLooping = ( thermalMinimum[ \"PotentialValue\" ]\n"
"                                < ( FunctionFromDictionary(\n"
"                                                   self.EffectivePotential,\n"
"                                                            fieldOrigin,\n"
"                                              evaporationTemperatureGuess )\n"
"                                - self.potentialScaler.functionAtOrigin ) )\n"
"            else:\n"
"                thermalArray = FieldDictionaryToArray( thermalMinimum[\n"
"                                                        \"FieldValues\" ] )\n"
"                keepLooping = ( numpy.sum( thermalArray**2 )\n"
"                                > self.minimumSeparationSquared )\n"
"# This relies on the field origin being a local minimum at the evaporation\n"
"# temperature - it should never be that MINUIT rolls the minimum to near\n"
"# the origin without the exact origin still being deeper.\n"
"        justBelowEvaporationTemperature = ( 0.5\n"
"                                            * evaporationTemperatureGuess )\n"
"        justAboveEvaporationTemperature = evaporationTemperatureGuess\n"
"        for narrowingStep in range( logTwoAccuracy ):\n"
"            evaporationTemperatureGuess = math.sqrt(\n"
"                                            justBelowEvaporationTemperature\n"
"                                        * justAboveEvaporationTemperature )\n"
"            print( \"Trying \"\n"
"                   + str( evaporationTemperatureGuess )\n"
"                   + \" GeV.\\n\" )\n"
"            tooLowTemperature = True\n"
"            self.SetTemperature( evaporationTemperatureGuess )\n"
"            thermalMinimum = self.TryToMinimize( minimumAtZeroTemperature[\n"
"                                                        \"FieldValues\" ] )\n"
"            if criticalComparedToOriginRatherThanEvaporation:\n"
"                tooLowTemperature = ( thermalMinimum[ \"PotentialValue\" ]\n"
"                                       < ( FunctionFromDictionary(\n"
"                                                   self.EffectivePotential,\n"
"                                                               fieldOrigin,\n"
"                                              evaporationTemperatureGuess )\n"
"                                - self.potentialScaler.functionAtOrigin ) )\n"
"            else:\n"
"                thermalArray = FieldDictionaryToArray( thermalMinimum[\n"
"                                                        \"FieldValues\" ] )\n"
"                tooLowTemperature = ( numpy.sum( thermalArray**2 )\n"
"                                      > self.minimumSeparationSquared )\n"
"            if tooLowTemperature:\n"
"                justBelowEvaporationTemperature = (\n"
"                                              evaporationTemperatureGuess )\n"
"            else:\n"
"                justAboveEvaporationTemperature = (\n"
"                                              evaporationTemperatureGuess )\n"
"        print( \"Evaporation temperature lies between \"\n"
"               + str( justBelowEvaporationTemperature )\n"
"               + \" GeV and \"\n"
"               + str( justAboveEvaporationTemperature )\n"
"               + \" GeV.\" )\n"
"        return [ justBelowEvaporationTemperature,\n"
"                 justAboveEvaporationTemperature ]\n"
"\n"
"\n"
"    def FindOptimalTunnelingTemperature( self,\n"
"                                         temperaturesForFit,\n"
"                                         justBelowCriticalTemperature,\n"
"                                         justAboveCriticalTemperature,\n"
"                                         dsbEvaporationTemperature ):\n"
"        guessedCriticalTemperature = math.sqrt(\n"
"                                               justBelowCriticalTemperature\n"
"                                           * justAboveCriticalTemperature )\n"
"        temperaturesWithThermalActions = []\n"
"        for trialTemperature in temperaturesForFit:\n"
"            optimalTunnelingTemperature = trialTemperature\n"
"            self.SetTemperature( trialTemperature )\n"
"            [ thermalDsbVacuum,\n"
"              thermalPanicVacuum ] = self.RollVacuaAtTemperature(\n"
"                                                          trialTemperature,\n"
"                                                dsbEvaporationTemperature )\n"
"            if not self.TunnelingPossible( thermalDsbVacuum,\n"
"                                           thermalPanicVacuum ):\n"
"                break\n"
"            thermalAction = self.CalculateAction(\n"
"                                            falseVacuum = thermalDsbVacuum,\n"
"                                           trueVacuum = thermalPanicVacuum,\n"
"                                                        deformPath = False,\n"
"                                                 thermalNotQuantum = True )\n"
"            alreadyExcluded = ( ( ( thermalAction / trialTemperature )\n"
"                                  + math.log( thermalAction ) )\n"
"                            < self.thermalActionOverTemperatureComparison )\n"
"            if alreadyExcluded:\n"
"# If the action at a given temperature is low enough that the DSB vacuum\n"
"# (or the symmetric vacuum that becomes the DSB vacuum through a 2nd-order\n"
"# phase transition) is unlikely to have survived the time when the Universe\n"
"# was at about that temperature, no further calculation is made.\n"
"                print( \"Parameter point is excluded by thermal tunneling\"\n"
"                       + \" at temperature \"\n"
"                       + str( trialTemperature )\n"
"                       + \" GeV.\" )\n"
"                return [ trialTemperature,\n"
"                         True,\n"
"                         str( self.ThermalTunnelingSurvivalProbability(\n"
"                                                          trialTemperature,\n"
"                                                        thermalAction ) ) ]\n"
"            print( \"\\nDirect path 3-dimensional action at temperature \"\n"
"                   + str( trialTemperature )\n"
"                   + \" GeV is \"\n"
"                   + str( thermalAction )\n"
"                   + \" GeV.\\n\\n\" )\n"
"            temperaturesWithThermalActions.append( [ trialTemperature,\n"
"                                                     thermalAction ] )\n"
"        if ( len( temperaturesWithThermalActions ) is 0 ):\n"
"            print( \"DSB vacuum seems to be lower than panic vacuum for\"\n"
"                   + \" each temperature tried: \"\n"
"                   + str( temperaturesForFit )\n"
"                   + \" (in GeV). No further thermal tunneling will be\"\n"
"                   + \"attempted.\" )\n"
"            self.ShouldTunnelThermally = False\n"
"            return [ 1.0,\n"
"                     False,\n"
"                     \"1.0\" ]\n"
"        print( \"temperaturesWithThermalActions = \"\n"
"               + str( temperaturesWithThermalActions ) )\n"
"        actionArray = numpy.array( [ temperatureWithThermalAction[ 1 ]\n"
"                                     for temperatureWithThermalAction in\n"
"                                     temperaturesWithThermalActions ] )\n"
"        transformationMatrix = numpy.array(\n"
"                [ numpy.array( [ ( ( trialTemperature**degreeOfPolynomial )\n"
"                                   / ( guessedCriticalTemperature\n"
"                                       - trialTemperature )**2 )\n"
"           for degreeOfPolynomial in range( len( temperaturesForFit ) ) ] )\n"
"                             for trialTemperature in temperaturesForFit ] )\n"
"        coefficientArray = numpy.linalg.solve( transformationMatrix,\n"
"                                               actionArray )\n"
"        print( \"3-dimensional action parameterized as ( \"\n"
"               + str( guessedCriticalTemperature )\n"
"               + \" - T )^( -2 ) * ( 0.0\" )\n"
"        for coefficientIndex in range( len( coefficientArray ) ):\n"
"            print( \"+ ( \"\n"
"                   + str( coefficientArray[ coefficientIndex ] )\n"
"                   + \" * T^\"\n"
"                   + str( coefficientIndex )\n"
"                   + \" )\" )\n"
"        print( \" ) GeV.\\n\\n\" )\n"
"        def ActionOverTMinusLogTFit( T ):\n"
"            returnValue = 0.0\n"
"            for coefficientIndex in range( len( coefficientArray ) ):\n"
"                returnValue += ( coefficientArray[ coefficientIndex ]\n"
"                                 * T**coefficientIndex )\n"
"            return ( ( returnValue\n"
"                       / ( ( guessedCriticalTemperature - T )**2 * T ) )\n"
"                     - 4.0 * math.log( T ) )\n"
"        temperatureMinimizer = minuit.Minuit( ActionOverTMinusLogTFit )\n"
"        temperatureMinimizer.limits[ \'T\' ] = ( 1.0,\n"
"                                             justBelowCriticalTemperature )\n"
"        optimalTunnelingTemperature = ( 0.5\n"
"                                        * justBelowCriticalTemperature )\n"
"        temperatureMinimizer.values[ \'T\' ] = optimalTunnelingTemperature\n"
"        try:\n"
"            temperatureMinimizer.migrad()\n"
"            optimalTunnelingTemperature = temperatureMinimizer.values[\n"
"                                                                    \'T\' ]\n"
"            print( \"Estimate based on -T^( -4 ) exp( S_3( T ) / T ) of\"\n"
"                   + \" temperature with lowest survival probability is \"\n"
"                   + str( optimalTunnelingTemperature )\n"
"                   + \" GeV.\\n\\n\" )\n"
"        except minuit.MinuitError as minuitError:\n"
"            print( \"PyMinuit threw an exception when trying to minimize\"\n"
"                   + \" the thermal survival probability. Half the\"\n"
"                   + \" critical temperature will now be used.\" )\n"
"            optimalTunnelingTemperature = ( 0.5\n"
"                                           * justBelowCriticalTemperature )\n"
"        except Exception as unexpectedException:\n"
"            print( \"Some lower-level exception happened:\\n\"\n"
"                   + str( unexpectedException )\n"
"                   + \"\\nwhile PyMinuit was trying to minimize the\"\n"
"                   + \" thermal survival probability. Half the critical\"\n"
"                   + \" temperature will now be used.\" )\n"
"            optimalTunnelingTemperature = ( 0.5\n"
"                                           * justBelowCriticalTemperature )\n"
"        self.SetTemperature( optimalTunnelingTemperature )\n"
"        [ thermalDsbVacuum,\n"
"          thermalPanicVacuum ] = self.RollVacuaAtTemperature(\n"
"                                               optimalTunnelingTemperature,\n"
"                                                dsbEvaporationTemperature )\n"
"        thermalAction = self.CalculateAction(\n"
"                                            falseVacuum = thermalDsbVacuum,\n"
"                                           trueVacuum = thermalPanicVacuum,\n"
"                                              deformPath = False,\n"
"                                              thermalNotQuantum = True )\n"
"        print( \"temperatureMinimizer.fval = \"\n"
"               + str( temperatureMinimizer.fval ) )\n"
"        print( \"Direct path tunneling at temperature \"\n"
"               + str( optimalTunnelingTemperature )\n"
"               + \" GeV:\\nEstimated 3-dimensional action = \"\n"
"               +str( optimalTunnelingTemperature\n"
"                     * ( temperatureMinimizer.fval\n"
"                         + ( 4.0\n"
"                            * math.log( optimalTunnelingTemperature ) ) ) )\n"
"               + \" GeV, actual = \"\n"
"               + str( thermalAction )\n"
"               + \" GeV.\" )\n"
"        return [ optimalTunnelingTemperature,\n"
"                 ( ( ( thermalAction / optimalTunnelingTemperature )\n"
"                     + math.log( thermalAction ) )\n"
"                   < self.thermalActionOverTemperatureComparison ),\n"
"                 str( self.ThermalTunnelingSurvivalProbability(\n"
"                                               optimalTunnelingTemperature,\n"
"                                                        thermalAction ) ) ]\n"
"\n"
"\n"
"    def ThermalTunnelingSurvivalProbability( self,\n"
"                                             exclusionTemperature,\n"
"                                             currentAction ):\n"
"        thermalExponent = ( lnOfThermalIntegrationFactor\n"
"                            - ( currentAction / exclusionTemperature ) )\n"
"        if ( thermalExponent < -exponentCutOff ):\n"
"            print( \"Warning! Thermal tunneling survival probability for\"\n"
"                   \" a very large ratio (\"\n"
"                   + str( currentAction / exclusionTemperature )\n"
"                   + \") of 3-dimensional action (\"\n"
"                   + str( currentAction )\n"
"                   + \") to temperature (\"\n"
"                   + str( exclusionTemperature )\n"
"                   + \") asked for, which would probably cause an\"\n"
"                   + \" overflow in the evaluation of the exponent.\"\n"
"                   + \" Returning a survival probability of 1.\" )\n"
"            return 1.0\n"
"# We don't need to consider thermalExponent > exponentCutOff as\n"
"# currentAction and exclusionTemperature should both be positive, and\n"
"# lnOfThermalIntegrationFactor < exponentCutOff.\n"
"        else:\n"
"            integratedDecayWidth = ( math.exp( thermalExponent )\n"
"                                     / currentAction )\n"
"            if ( integratedDecayWidth > exponentCutOff ):\n"
"                print( \"Warning! Integrated thermal decay width is so\"\n"
"                       + \" large and positive that an overflow error\"\n"
"                       + \" would occur when exponentiating to get the\"\n"
"                       + \" survival probability, so returning a survival\"\n"
"                       + \" probability of 0.\" )\n"
"                return 0.0\n"
"# We don't need to consider integratedDecayWidth < -exponentCutOff as\n"
"# math.exp( thermalExponent ) must be positive, and currentAction should be\n"
"# positive.\n"
"            return math.exp( -integratedDecayWidth )\n"
"\n"
"# End of Vevacious class\n"
"\n";
    outputFile.close();
  }

  void VevaciousRunner::appendResultsToSlha( std::string const& slhaFilename,
                          bool const properFormatRatherThanSspReadable ) const
  // this takes the XML-format results from the file called resultsFilename,
  // & appends the results & warnings in custom SLHA blocks to the end of the
  // file name slhaFilename.
  {
    double stabiltyResult( -99.0 );
    // -99 will do as an error code for the SLHA block for the moment.
    std::string stabilityVerdict( "error" );
    std::map< std::string, std::string > inputMinimumDepthAndVevs;
    std::map< std::string, std::string > globalMinimumDepthAndVevs;
    double calculatedLifetime( -99.0 );
    std::string thermalSurvival( "-99.0" );
    double exclusionTemperature( -99.0 );
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
          std::map< std::string, std::string >::const_iterator
          slhaFinder( resultParser.getCurrentElementAttributes().find(
                                                       "slha_summary_code" ) );
          if( resultParser.getCurrentElementAttributes().end()
              != slhaFinder )
          {
            stabiltyResult = BOL::StringParser::stringToDouble(
                                      BOL::StringParser::trimFromFrontAndBack(
                                                            slhaFinder->second,
                                                                    "\'\"" ) );
          }
        }
        else if( resultParser.currentElementNameMatches(
                                                        "quantum_stability" ) )
        {
          stabilityVerdict.assign(
                              resultParser.getTrimmedCurrentElementContent() );
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
          calculatedLifetime = doubleFromQuotedString(
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
        else if( resultParser.currentElementNameMatches(
                                                    "tunneling_temperature" ) )
        {
          exclusionTemperature = doubleFromQuotedString(
                              resultParser.getTrimmedCurrentElementContent() );
          std::map< std::string, std::string >::const_iterator
          actionFinder( resultParser.getCurrentElementAttributes().find(
                                                      "survival_probability" ) );
          if( resultParser.getCurrentElementAttributes().end()
              != actionFinder )
          {
            thermalSurvival
            = BOL::StringParser::trimFromFrontAndBack( actionFinder->second,
                                                       "\'\"" );
          }
        }
        else if( resultParser.currentElementNameMatches( "warning" ) )
        {
          warningLine.assign( resultParser.getTrimmedCurrentElementContent() );
          Vevacious::SarahInterpreter::removeNewlinesFrom( warningLine );
          warningBlock.append( " " );
          warningBlock.append( slhaIndexMaker.intToString( ++warningNumber ) );
          warningBlock.append( "    " );
          warningBlock.append( warningLine );
          warningBlock.append( "\n" );
        }
      }
      resultParser.closeFile();

      if( stabiltyResult == -2.0 )
      {
        stabilityVerdict.assign( "long-lived_but_thermally_excluded" );
      }
      else if( stabiltyResult == -3.0 )
      {
        stabilityVerdict.assign( "long-lived_but_maybe_thermally_excluded" );
      }

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
      << "BLOCK VEVACIOUSRESULTS # results from Vevacious version "
      << vevaciousVersion << ", documented in " << vevaciousDocumentation
      << "\n"
      << writeVevaciousResultsBlockLine( 0,
                                         0,
                                         stabiltyResult,
                                         stabilityVerdict,
                                         properFormatRatherThanSspReadable )
      << "    # stability of input\n"
      << writeVevaciousResultsBlockLine( 0,
                                         1,
                                         calculatedLifetime,
                                         actionCalculation,
                                         properFormatRatherThanSspReadable )
      << "    # tunneling time in Universe ages / calculation type\n"
      << writeVevaciousResultsBlockLine( 0,
                                         2,
                                         exclusionTemperature,
                                         thermalSurvival,
                                         properFormatRatherThanSspReadable )
      << "    # estimated best tunneling temperature"
      << " / survival probability at this temperature\n";
      int vevNumber( -1 );
      for( std::map< std::string, std::string >::const_iterator
           whichVev( inputMinimumDepthAndVevs.begin() );
           inputMinimumDepthAndVevs.end() != whichVev;
           ++whichVev )
      {
        outputFile
        << writeVevaciousResultsBlockLine( 1,
                                           (++vevNumber),
                                    doubleFromQuotedString( whichVev->second ),
                                           whichVev->first,
                                           properFormatRatherThanSspReadable );
        if( inputMinimumDepthAndVevs.begin() == whichVev )
        {
          outputFile
          << "    # DSB vacuum potential energy\n";
        }
        else
        {
          outputFile
          << "    # DSB vacuum VEV\n";
        }
      }
      vevNumber = -1;
      for( std::map< std::string, std::string >::const_iterator
           whichVev( globalMinimumDepthAndVevs.begin() );
           globalMinimumDepthAndVevs.end() != whichVev;
           ++whichVev )
      {
        outputFile
        << writeVevaciousResultsBlockLine( 2,
                                           (++vevNumber),
                                    doubleFromQuotedString( whichVev->second ),
                                           whichVev->first,
                                           properFormatRatherThanSspReadable );
        if( globalMinimumDepthAndVevs.begin() == whichVev )
        {
          outputFile
          << "    # panic vacuum potential\n";
        }
        else
        {
          outputFile
          << "    # panic vacuum VEV\n";
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
        << (*whichSet)[ whichField ] << " )";
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
"# (B.O'L. would like to apologize about the state of this file: ideally it\n"
"# would be a lot neater, and respect modern programming practices. Maybe a\n"
"# later version will tidy it up properly...)\n"
"#\n"
"#\n"
"from __future__ import division\n"
"import math\n"
"import numpy\n"
"import minuit\n"
"import VevaciousParameterDependent as VPD\n"
"import VevaciousTreeLevelExtrema as VTE\n\n"
"\n"
"treeLevelExtrema = VTE.pointsToTry\n"
"vcs = VPD.Vevacious(\n"
"                EffectivePotential = VPD.LoopAndThermalCorrectedPotential )\n"
"# Various settings can be changed here: e.g. to use the tree-level\n"
"# potential for minimizing and tunneling, one could use\n"
"# EffectivePotential = VPD.TreeLevelPotential\n"
"# Even a custom potential function can be inserted here, as long as it\n"
"# takes the correct arguments. Have a look at\n"
"# VevaciousParameterDependent.py to see the correct format.\n"
"# IF YOU DO USE THE TREE-LEVEL POTENTIAL, PLEASE REMEMBER TO TURN OFF\n"
"# THERMAL CALCULATIONS! The thermal corrections are formally a loop effect,\n"
"# and the tree-level potential is independent of the temperature. This\n"
"# means that Vevacious would spend forever trying to find the critical\n"
"# temperature at which tunneling from the DSB minimum to the panic vacuum\n"
"# becomes impossible This can just be set by uncommenting the next line.\n"
"#vcs.shouldDoThermalCalculations = False\n"
"vcs.allowedRunningTime = 3600.0\n"
"# Allowing an hour of running time is maybe excessive...\n"
"vcs.SetFieldValueLimit( 5.0 * VPD.energyScale )\n"
"# MINUIT is now not allowed to roll outside a hypercube defined by all the\n"
"# field values being less than five times the renormalization scale."
"\n"
"# WriteExtrema is not important for deciding whether a parameter point is\n"
"# metastable, but it helps with debugging.\n"
"vcs.WriteExtrema( treeLevelExtrema, \"Vevacious_tree-level_extrema.txt\" )\n"
"\n"
"# RollExtrema includes the automatic nudging off saddle points according to\n"
"# vcs.nudgeList.\n"
"# It also sorts the minima and sets vcs.panicVacuum and vcs.globalMinimum,\n"
"# and vcs.dsbVacuumIsMetastable.\n"
"vcs.RollExtrema( treeLevelExtrema )\n"
"print( \"DSB vacuum: \" + vcs.ExtremumAsMathematica( vcs.dsbVacuum ) )\n"
"\n"
"# WriteMinima is also not important for deciding whether a parameter point\n"
"# is metastable, but it helps with debugging.\n"
"vcs.WriteMinima( \"Vevacious_loop-corrected_minima.txt\" )\n"
"\n"
"metastabilityVerdict = \"not_calculated\"\n"
"quantumTunnelingActionType = \"not_calculated\"\n"
"quantumStabilityVerdict = \"not_calculated\"\n"
"thermalStabilityVerdict = \"not_calculated\"\n"
"alreadyExcluded = False\n"
"exclusionTemperatureString = \"not_calculated\"\n"
"quantumTunnelingTimeInUniverseAgesString = \"not_calculated\"\n"
"thermalTunnelingSurvivalProbabilityString = \"not_calculated\"\n"
"slhaSummaryCode = -99\n"
"\n"
"\n"
"if ( not vcs.dsbVacuumIsMetastable ):\n"
"    print( \"DSB vacuum is stable (as far as the model file allows).\" )\n"
"    metastabilityVerdict = \"stable\"\n"
"    quantumTunnelingActionType = \"unnecessary\"\n"
"    quantumStabilityVerdict = \"stable\"\n"
"    thermalStabilityVerdict = \"stable\"\n"
"    exclusionTemperatureString = \"unnecessary\"\n"
"    quantumTunnelingTimeInUniverseAgesString = \"-1.0\"\n"
"    thermalTunnelingSurvivalProbabilityString = \"1.0\"\n"
"    slhaSummaryCode = 1\n"
"else:\n"
"# The deepest minimum found by vcs.RollExtrema is recorded as\n"
"# vcs.globalMinimum and it may happen to correspond to the DSB minimum.\n"
"# If there were minima found which are deeper than the DSB minimum, the one\n"
"# closest to the DSB minimum (or reflections of it in any field value, as\n"
"# the potential may be invariant under a set of sign flips of the fields)\n"
"# is recorded as the panic vacuum, in vcs.panicVacuum. The panic vacuum\n"
"# does not necessarily correspond to vcs.globalMinimum. Strictly, all\n"
"# minima should be checked for tunneling, but it is computationally\n"
"# expensive, and it seems that for points of phenomenological interest, the\n"
"# bubble wall is the dominant term so the closest minima are more likely to\n"
"# have the lowest actions.\n"
"    print( \"Panic vacuum: \"\n"
"           + vcs.ExtremumAsMathematica( vcs.panicVacuum ) )\n"
"    print( \"Global minimum found: \"\n"
"           + vcs.ExtremumAsMathematica( vcs.globalMinimum )\n"
"           + \"\\n\\n\" )\n"
"    metastabilityVerdict = \"metastable\"\n"
"    slhaSummaryCode = 0\n"
"\n"
"\n"
"if ( vcs.dsbVacuumIsMetastable\n"
"     and\n"
"     vcs.ShouldTunnel ):\n"
"# First, a direct path between the minima is taken to get an upper bound on\n"
"# the bounce action.\n"
"    quantumTunnelingActionType = \"direct_path\"\n"
"    quantumStabilityVerdict = \"long-lived\"\n"
"    thermalStabilityVerdict = \"high_survival_probability\"\n"
"    vcs.SetTemperature( 0.0 )\n"
"    quantumAction = vcs.CalculateAction( falseVacuum = vcs.dsbVacuum,\n"
"                                         trueVacuum = vcs.panicVacuum,\n"
"                                         deformPath = False,\n"
"                                         thermalNotQuantum = False )\n"
"    print( \"Direct path zero-temperature 4-dimensional action is \"\n"
"           + str( quantumAction )\n"
"           + \"\\n\\n\" )\n"
"    quantumTunnelingTimeInUniverseAgesString = str(\n"
"                      vcs.QuantumTunnelingTimeInInverseGev( quantumAction )\n"
"                                     / vcs.ageOfKnownUniverseInInverseGev )\n"
"    alreadyExcluded = ( quantumAction <= vcs.quantumActionThreshold )\n"
"    if alreadyExcluded:\n"
"        quantumStabilityVerdict = \"short-lived\"\n"
"        thermalStabilityVerdict = \"low_survival_probability\"\n"
"        slhaSummaryCode = -1\n"
"    elif vcs.ShouldTunnelThermally:\n"
"        optimalTunnelingTemperature = 0.0\n"
"# If the parameter point is not excluded by the naive straight path at zero\n"
"# temperature, then we check thermal tunneling by direct paths to see if\n"
"# a direct path can exclude the point, & to get an estimate of what\n"
"# temperature to use for a full calculation with path deformation.\n"
"        print( \"Upper bound on zero-temperature full effective action by\"\n"
"               + \" direct path is too high to exclude point, so trying\"\n"
"               + \" direct paths at non-zero temperatures.\" )\n"
"        dsbEvaporationTemperature = vcs.FindEvaporationTemperature(\n"
"                                                             vcs.dsbVacuum,\n"
"                                                                    2,\n"
"                                                               False )[ 1 ]\n"
"        if not ( dsbEvaporationTemperature > 0.0 ):\n"
"            thermalStabilityVerdict = \"low_survival_probability\"\n"
"            exclusionTemperatureString = str( dsbEvaporationTemperature )\n"
"            alreadyExcluded = True\n"
"            dsbDepth = VPD.FunctionFromDictionary( vcs.EffectivePotential,\n"
"                                          vcs.dsbVacuum[ \"FieldValues\" ],\n"
"                                                   0.0 )\n"
"            originDepth = VPD.FunctionFromDictionary(\n"
"                                                    vcs.EffectivePotential,\n"
"                                                      VPD.fieldOrigin,\n"
"                                                      0.0 )\n"
"            vcs.LogWarning( \"Point with all fields = 0 is deeper (\"\n"
"                            + str( originDepth )\n"
"                            + \" GeV^4) than DSB vacuum (\"\n"
"                            + str( dsbDepth )\n"
"                            + \" GeV^4) at T = 0!\" )\n"
"        else:\n"
"            criticalTemperatureRange = vcs.FindEvaporationTemperature(\n"
"                                                           vcs.panicVacuum,\n"
"                                                                       8,\n"
"                                                                     True )\n"
"# The critical temperature is above criticalTemperatureRange[ 0 ] and below\n"
"# criticalTemperatureRange[ 1 ], and the range should be no more than\n"
"# 2^(-logTwoAccuracy) times the critical temperature in extent.\n"
"            fitNodes = 3\n"
"            lowestFitTemperature = dsbEvaporationTemperature\n"
"# We assume that the optimal tunneling temperature will be above the DSB\n"
"# evaporation temperature, but if the DSB evaporation temperature is close\n"
"# enough to the critical temperature, we include temperatures below it for\n"
"# the fit, even though it introduces a kink which is actually pretty badly\n"
"# approximated by a polynomial. (If the DSB evaporation temperature is\n"
"# higher than the panic vacuum's evaporation temperature, then there is no\n"
"# kink and we consider a lower resolution polynomial, across the range\n"
"# from 1 GeV up to 0.9 times the panic evaporation temperature.)\n"
"            if ( ( criticalTemperatureRange[ 0 ]\n"
"                   / dsbEvaporationTemperature ) < 2.0 ):\n"
"                lowestFitTemperature = 1.0\n"
"                if ( criticalTemperatureRange[ 0 ]\n"
"                     > dsbEvaporationTemperature ):\n"
"                    fitNodes = 10\n"
"            temperatureFitDifference = ( 0.9\n"
"                                         * ( criticalTemperatureRange[ 0 ]\n"
"                                            - dsbEvaporationTemperature ) )\n"
"            temperaturesForFit = [ ( lowestFitTemperature\n"
"                                     + ( ( stepIndex\n"
"                                           * temperatureFitDifference )\n"
"                                         / ( fitNodes - 1.0 ) ) )\n"
"                                   for stepIndex in range( fitNodes ) ]\n"
"            [ optimalTunnelingTemperature,\n"
"              alreadyExcluded,\n"
"              thermalTunnelingSurvivalProbabilityString\n"
"                                   ] = vcs.FindOptimalTunnelingTemperature(\n"
"                                                        temperaturesForFit,\n"
"                                             criticalTemperatureRange[ 0 ],\n"
"                                             criticalTemperatureRange[ 1 ],\n"
"                                                dsbEvaporationTemperature )\n"
"            exclusionTemperatureString = str(\n"
"                                              optimalTunnelingTemperature )\n"
"# vcs.FindOptimalTunnelingTemperature fits n = len( temperaturesForFit )\n"
"# coefficents a_1 to a_n for a function\n"
"# a_1 (Tc - T )^(-2) + a_2 T^0 + a_3 T + ...\n"
"# to approximate the temperature (T) dependence (where Tc is the criticial\n"
"# temperature) of the 3-dimensional thermal action calculated by a direct\n"
"# path in field space, which is then minimized by PyMinuit to find an\n"
"# estimate of the optimal tunneling temperature. It also returns True if\n"
"# the cautiously-high estimate of the survival probability at any of the\n"
"# given temperatures is below the threshold given by the Vevacious\n"
"# initialization XML file.\n"
"        if alreadyExcluded:\n"
"            thermalStabilityVerdict = \"low_survival_probability\"\n"
"            slhaSummaryCode = -2\n"
"\n"
"\n"
"    if ( vcs.ShouldDeformTunnelPaths\n"
"         and\n"
"         ( not ( alreadyExcluded\n"
"                 or\n"
"                 vcs.AllowedRunningTimeExceeded() ) ) ):\n"
"# Here we continue to check the zero-temperature quantum tunneling (if the\n"
"# direct path at zero temperature did not already exclude the point):\n"
"        vcs.SetTemperature( 0.0 )\n"
"        quantumAction = vcs.CalculateAction( falseVacuum = vcs.dsbVacuum,\n"
"                                             trueVacuum = vcs.panicVacuum,\n"
"                                             deformPath = True,\n"
"                                             thermalNotQuantum = False,\n"
"                                             innerLoopMaxDeformations = 10,\n"
"                                            outerLoopMaxDeformations = 10 )\n"
"# CosmoTransitions nests a loop of deformations inside a loop for some\n"
"# unclear reason. The arguments passed here are intentionally low, as the\n"
"# path deformation takes a _very_ long time, and the returns on further\n"
"# iterations usually diminish very quickly.\n"
"        quantumTunnelingActionType = \"deformed_path\"\n"
"        print( \"Final deformed path zero-temperature 4-dimensional\"\n"
"               + \" action is \"\n"
"               + str( quantumAction )\n"
"               + \" units of h-bar.\\n\\n\" )\n"
"        quantumTunnelingTimeInInverseGev = (\n"
"                                      vcs.QuantumTunnelingTimeInInverseGev(\n"
"                                                          quantumAction ) )\n"
"        print( \"Zero-temperature tunneling time in 1/GeV = \"\n"
"               + str( quantumTunnelingTimeInInverseGev ) )\n"
"        print( \"Age of known Universe in 1/GeV = \"\n"
"               + str( vcs.ageOfKnownUniverseInInverseGev ) )\n"
"        quantumTunnelingTimeInUniverseAges = (\n"
"                                           quantumTunnelingTimeInInverseGev\n"
"                                     / vcs.ageOfKnownUniverseInInverseGev )\n"
"        quantumTunnelingTimeInUniverseAgesString = str(\n"
"                                       quantumTunnelingTimeInUniverseAges )\n"
"        print( \"Tunneling time in Universe-ages = \"\n"
"               + quantumTunnelingTimeInUniverseAgesString )\n"
"        quantumTunnelingTimeInSeconds = str(\n"
"                                            vcs.ageOfKnownUniverseInSeconds\n"
"                                     * quantumTunnelingTimeInUniverseAges )\n"
"        print( \"Zero-temperature tunneling time estimate is \"\n"
"               + quantumTunnelingTimeInSeconds\n"
"               + \" seconds (age of known Universe = \"\n"
"               + str( vcs.ageOfKnownUniverseInSeconds )\n"
"               + \" seconds).\\n\\n\" )\n"
"        if ( quantumAction < vcs.quantumActionThreshold ):\n"
"            quantumStabilityVerdict = \"short-lived\"\n"
"            slhaSummaryCode = -1\n"
"        elif ( vcs.ShouldTunnelThermally\n"
"               and\n"
"               ( not vcs.AllowedRunningTimeExceeded() ) ):\n"
"# If the parameter point was not excluded just by quantum tunneling, we\n"
"# calculate whether fully deformed thermal tunneling excludes it.\n\n"
"            vcs.SetTemperature( optimalTunnelingTemperature )\n"
"            print( \"Deforming tunneling path at temperature \"\n"
"                   + str( optimalTunnelingTemperature )\n"
"                   + \" GeV.\" )\n"
"            print( \"At this temperature, the threshold 3-dimensional\"\n"
"                   + \" action is \"\n"
"                   + str( optimalTunnelingTemperature\n"
"                          * ( vcs.thermalActionOverTemperatureComparison\n"
"                              - math.log( optimalTunnelingTemperature ) ) )\n"
"                   + \" GeV (aggressive; a cautious threshold is less by\"\n"
"                   + \" roughly ten times \"\n"
"                   + str( optimalTunnelingTemperature )\n"
"                   + \" GeV).\\n\\n\" )\n"
"            thermalAction = vcs.CalculateAction(\n"
"                                               falseVacuum = vcs.dsbVacuum,\n"
"                                              trueVacuum = vcs.panicVacuum,\n"
"                                                deformPath = True,\n"
"                                                thermalNotQuantum = True,\n"
"                                             innerLoopMaxDeformations = 10,\n"
"                                            outerLoopMaxDeformations = 10 )\n"
"            print( \"Final 3-dimensional action  at this temperature = \"\n"
"                   + str( thermalAction )\n"
"                   + \" GeV.\\n\\n\" )\n"
"            if not ( thermalAction > 0.0 ):\n"
"                thermalStabilityVerdict = \"unsure_survival_probability\"\n"
"                slhaSummaryCode = -3\n"
"                print( \"CosmoTransitions seems to have failed to\"\n"
"                       + \" converge: marking the survival probability as\"\n"
"                       + \"uncertain.\\n\\n\" )\n"
"            if ( ( ( thermalAction / optimalTunnelingTemperature )\n"
"                     + math.log( thermalAction ) )\n"
"                   < vcs.thermalActionOverTemperatureComparison ):\n"
"                thermalStabilityVerdict = \"low_survival_probability\"\n"
"                slhaSummaryCode = -2\n"
"            elif ( ( ( thermalAction / optimalTunnelingTemperature )\n"
"                     + math.log( optimalTunnelingTemperature ) )\n"
"                   < vcs.thermalActionOverTemperatureComparison ):\n"
"                thermalStabilityVerdict = \"unsure_survival_probability\"\n"
"                slhaSummaryCode = -3\n"
"            print( \"Temperature for basing thermal tunneling exclusion\"\n"
"                   + \" = \"\n"
"                   + exclusionTemperatureString\n"
"                   + \" GeV.\\n\\n\" )\n"
"            thermalTunnelingSurvivalProbabilityString = str(\n"
"                                   vcs.ThermalTunnelingSurvivalProbability(\n"
"                                               optimalTunnelingTemperature,\n"
"                                                          thermalAction ) )\n"
"\n"
"\n"
"# Finally the output file is written:\n"
"outputText = ( \"  <reference version=\\\"" << vevaciousVersion << "\\\"\"\n"
"               + \" citation=\\\"" << vevaciousDocumentation
<<                                                             "\\\" />\\n\"\n"
"               + \"  <stability slha_summary_code=\\\"\"\n"
"               + str( slhaSummaryCode )\n"
"               + \"\\\" > \"\n"
"               + metastabilityVerdict\n"
"               + \" </stability>\\n\"\n"
"               + \"  <quantum_stability> \"\n"
"               + quantumStabilityVerdict\n"
"               + \" </quantum_stability>\\n\"\n"
"               + \"  <thermal_stability> \"\n"
"               + thermalStabilityVerdict\n"
"               + \" </thermal_stability>\\n\"\n"
"               + \"  <global_minimum   relative_depth=\\\"\"\n"
"               + str( vcs.panicVacuum[ \"PotentialValue\" ] )\n"
"               + \"\\\" \"\n"
"               + VPD.UserFieldsAsXml( vcs.panicVacuum[ \"FieldValues\" ] )\n"
"               + \" />\\n  <input_minimum   relative_depth=\\\"\"\n"
"               + str( vcs.dsbVacuum[ \"PotentialValue\" ] )\n"
"               + \"\\\" \"\n"
"               + VPD.UserFieldsAsXml( vcs.dsbVacuum[ \"FieldValues\" ] )\n"
"               + \" />\\n  <lifetime  action_calculation=\\\"\"\n"
"               + quantumTunnelingActionType\n"
"               + \"\\\" > \"\n"
"               + quantumTunnelingTimeInUniverseAgesString\n"
"               + \" </lifetime>\\n\"\n"
"               + \" <tunneling_temperature survival_probability=\\\"\"\n"
"               + thermalTunnelingSurvivalProbabilityString\n"
"               + \"\\\" > \"\n"
"               + exclusionTemperatureString\n"
"               + \" </tunneling_temperature>\" )\n"
"outputFile = open( VPD.outputFile, \"w\" )\n"
"outputFile.write( \"<Vevacious_result>\\n\"\n"
"                  + outputText )"
"# Each warning is printed as an XML element:\n"
"for warningMessage in vcs.warningMessages:\n"
"    outputFile.write( \"\\n  <warning>\\n  \"\n"
"                      + warningMessage\n"
"                      + \"\\n  </warning>\" )\n"
"outputFile.write( \"\\n</Vevacious_result>\\n\" )\n"
"outputFile.close()\n"
"print( \"Result summary (not recapping warnings):\\n\""
"       + outputText )\n"
"\n";
    outputFile.close();
  }

} /* namespace Vevacious */
