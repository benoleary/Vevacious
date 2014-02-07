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
  std::string const VevaciousRunner::vevaciousVersion( "1.1.00" );
  std::string const
  VevaciousRunner::vevaciousDocumentation( "arXiv:1307.1477 (hep-ph)" );
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
    lifetimeThreshold( 0.1 )
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
                                                               "0.1" ) ) )
  {
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
"reducedPlanckMass = 2.43E+18\n"
"numberOfHorizons = (1.0E+14)**3\n"
"lnOfNumberOfHorizons = math.log( numberOfHorizons )\n"
"energyScale = " << sarahInterpreter.getSlhaScale() << "\n"
"energyScaleFourth = ( energyScale**4 )\n"
"inverseScale = ( 1.0 / energyScale )\n"
"inverseScaleSquared = ( 1.0 / ( energyScale * energyScale ) )\n"
"inverseScaleFourthed = ( 1.0 / energyScaleFourth )\n"
"numericalStepSize = 1.0\n"
"exponentCutOff = math.log( 0.25 * sys.float_info.max )\n"
"\n"
"# Unfortunately, due to UTF32 issues, PyMinuit is restricted to\n"
"\n"
"# First some general functions:\n"
"\n"
"pathToCosmotransitions = \"" << pathToCosmotransitions << "\"\n"
"sys.path.append( pathToCosmotransitions )\n"
"import pathDeformation as CTPD\n"
"import finiteT as CTFT\n"
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
"        if( ( 1.0E-6 * inverseScaleSquared ) < massSquaredMagnitude ):\n"
"            summedCorrection += ( massSquaredMagnitude\n"
"                                  * massSquaredMagnitude \n"
"                                  * ( math.log( massSquaredMagnitude\n"
"                                                * inverseScaleSquared )\n"
"                                      - subtractionConstant ) )\n"
"    return ( overallFactor * summedCorrection )\n"
"\n"
"loopFactor = ( 1.0 / ( 64.0 * math.pi * math.pi ) )\n"
"\n"
"\n"
"# The following are approximations to the J_{+} & J_{-} functions:\n"
"# J_{+}( r ) = integral from 0 to infinity by dx of \n"
"# x^2 ln( 1 - exp( -sqrt( x^2 + r^2 ) ) )\n"
"# J_{-}( r ) = integral from 0 to infinity by dx of \n"
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
"\n"
"        self.EffectivePotential = EffectivePotential\n"
"        self.potentialScaler = PotentialScaler( self.EffectivePotential,\n"
"                                    temperatureValue = currentTemperature )\n"
"        self.PyMinuit = minuit.Minuit(\n"
"                   self.potentialScaler.ScaledFunctionFromScaledArguments )\n"
"        self.rollingToleranceSquared = (0.1)**2\n"
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
"        self.dsbVacuum = self.TryToMinimize( dsbInput )\n"
"        self.dsbArray = FieldDictionaryToArray( self.dsbVacuum[\n"
"                                                        \"FieldValues\" ] )\n"
"        self.panicVacuum = self.dsbVacuum\n"
"        self.panicArray = self.dsbArray\n"
"        self.globalMinimum = self.dsbVacuum\n"
"        self.globalArray = self.dsbArray\n"
"\n"
"        self.ageOfKnownUniverseInSeconds = ageOfKnownUniverseInSeconds\n"
"        self.ageOfKnownUniverseInInverseGev = ( ageOfKnownUniverseInSeconds\n"
"                                                / 6.582119E-16 )\n"
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
"                                   * self.ageOfTheKnownUniverseInInverseGev\n"
"                                         * fourthRootOfSolitonicFactorA ) )\n"
"        self.thermalActionThreshold = 1000.0\n"
"        self.warningMessages = []\n"
"        self.startingTime = time.clock()\n"
"        self.allowedRunningTime = 100.0\n"
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
"            LogWarning( \"PyMinuit had problems starting at \"\n"
"                        + UserFieldsAsMathematica( fieldValueDictionary )\n"
"                        + \"! [minuit.MinuitError: \"\n"
"                        + str( minuitError )\n"
"                        + \"]. PyMinuit stopped at \"\n"
"                        + UserFieldsAsMathematica(\n"
"                        self.FieldValuesFromMinuit( minuitObject.values ) )\n"
"                        + \".  Minuit\'s estimate of how much deeper it\"\n"
"                        + \" should go is \"\n"
"                 + str( self.PotentialValueFromMinuit( minuitObject.edm ) )\n"
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
HERE!
"\n"
"    def TryToMinimizeIncludingRescaling( self, fieldValueDictionary ):\n"
"        foundExtremum = self.TryToMinimize( fieldValueDictionary )\n"
"        if ( self.CloseEnoughToDsbOrSignFlip( foundExtremum[ \"FieldValues\" ] )\n"
"             and\n"
"             not self.CloseEnoughToDsbOrSignFlip( fieldValueDictionary ) ):\n"
"            scaledValues = fieldValueDictionary.copy()\n"
"            for fieldKey in fieldValueDictionary.keys():\n"
"                scaledValues[ fieldKey ] *= self.fieldScaling\n"
"            scaledPoint = self.FieldDictionaryToExtremum( scaledValues )\n"
"            self.LogWarning( \"Starting point differed from DSB vacuum but rolled there.\"\n"
"                             + \" Trying again with scaled field configuration (\"\n"
"                             + self.ExtremumAsMathematica( self.FieldDictionaryToExtremum( scaledPoint[ \"FieldValues\" ] ) ) + \").\" )\n"
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
"            foundExtrema.append( self.TryToMinimizeIncludingRescaling( startingFieldConfiguration ) )\n"
"        if ( len( self.saddleSplitNudges ) > 0 ):\n"
"            for foundExtremum in foundExtrema:\n"
"                isMinimum, nudgeVector = self.CheckHessian( foundExtremum )\n"
"                if ( isMinimum ):\n"
"                    self.foundMinima.append( foundExtremum )\n"
"                else:\n"
"                    foundSaddlePoints.append( [ foundExtremum, nudgeVector ] )\n"
"        else:\n"
"            self.foundMinima = foundExtrema\n"
"# At this point, the field configurations from startingFieldConfigurations have\n"
"# been divided into minima in self.foundMinima and saddle points (& possibly\n"
"# maxima) in foundSaddlePoints, if self.saddleSplitNudges is not an empty list.\n"
"# (If self.saddleSplitNudges is an empty list, at this point self.foundMinima\n"
"# contains all the extrema from the 1st pass, & this function won\'t do anything\n"
"# else.)\n"
"        for saddleSplitNudge in self.saddleSplitNudges:\n"
"# For each requested nudging, nudgedFieldConfigurations is made from all the\n"
"# pairs of nearby configurations to each in foundSaddlePoints,\n"
"# foundSaddlePoints is then emptied, & each configuration is rolled & the\n"
"# rolled configuration is put into either self.foundMinima if a minimum or\n"
"# foundSaddlePoints again if not, & then the loop iterates.\n"
"            nudgedFieldConfigurations = []\n"
"            for foundSaddlePoint in foundSaddlePoints:\n"
"                nudgedFieldConfigurations.extend( self.NudgeOff( foundSaddlePoint, saddleSplitNudge ) )\n"
"            foundSaddlePoints = []\n"
"            for nudgedFieldConfiguration in nudgedFieldConfigurations:\n"
"                foundExtremum = self.TryToMinimizeIncludingRescaling( nudgedFieldConfiguration )\n"
"                isMinimum, nudgeVector = self.CheckHessian( foundExtremum )\n"
"                if ( isMinimum ):\n"
"                    self.foundMinima.append( foundExtremum )\n"
"                else:\n"
"                    foundSaddlePoints.append( [ foundExtremum, nudgeVector ] )\n"
"# At this point, any remaining saddle points / maxima are just appended to\n"
"# self.foundMinima so that they can still be tunneled to, regardless of the\n"
"# fact that the Universe would still evolve away from such configurations.\n"
"        for foundSaddlePoint in foundSaddlePoints:\n"
"            self.foundMinima.append( foundSaddlePoint[ 0 ] )\n"
"# Now the minima are sorted by depth:\n"
"        tunnelDistance = -1.0\n"
"        for foundMinimum in self.foundMinima:\n"
"            if ( foundMinimum[ \"PotentialValue\" ] < self.dsbVacuum[ \"PotentialValue\" ] ):\n"
"                self.dsbVacuumIsMetastable = True\n"
"                fieldArray = FieldDictionaryToArray( foundMinimum[ \"FieldValues\" ] )\n"
"                currentDistanceToDsb = numpy.sum( ( fieldArray - self.dsbArray )**2 )\n"
"                if ( ( tunnelDistance < 0.0 )\n"
"                     or\n"
"                     ( currentDistanceToDsb < tunnelDistance ) ):\n"
"                    tunnelDistance = currentDistanceToDsb\n"
"                    self.panicVacuum = foundMinimum\n"
"                    self.panicArray = fieldArray\n"
"                if ( foundMinimum[ \"PotentialValue\" ] < self.globalMinimum[ \"PotentialValue\" ] ):\n"
"                    self.globalMinimum = foundMinimum\n"
"                    self.globalArray = fieldArray\n"
"\n"
"\n"
"    def CheckHessian( self, foundExtremum ):\n"
"        isMinimum = True\n"
"        nudgeVector = list( originArray )\n"
"        eigensystemOfHessian = numpy.linalg.eigh( NumericalHessian( self.EffectivePotential,\n"
"                                                                    foundExtremum[ \"FieldValues\" ] ) )\n"
"# It doesn\'t matter what value mostNegativeEigenvalueValue has if it\'s positive\n"
"# because it is only there to record the lowest eigenvalue if it is less than\n"
"# or equal to zero. \n"
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
"            nudgeVector = eigensystemOfHessian[ 1 ][ :, mostNegativeEigenvalueIndex ]\n"
"        return [ isMinimum, nudgeVector ]\n"
"\n"
"\n"
"    def WriteMinima( self, fileToWrite ):\n"
"        outputFile = open( fileToWrite, \"w\" )\n"
"        outputFile.write( \"{\\n\" )\n"
"        for foundMinimum in self.foundMinima:\n"
"            outputFile.write( self.ExtremumAsMathematica( foundMinimum ) + \"\\n\" )\n"
"        outputFile.write( \"}\" )\n"
"        outputFile.close()\n"
"\n"
"\n"
"    def CalculateAction( self, \n"
"                         falseVacuum = None,\n"
"                         trueVacuum = None,\n"
"                         deformPath = True,\n"
"                         thermalNotQuantum = None,\n"
"                         **keywordArguments ):\n"
"        if ( falseVacuum is None ):\n"
"            falseVacuum = self.dsbVacuum\n"
"        if ( trueVacuum is None ):\n"
"            trueVacuum = self.panicVacuum\n"
"        if thermalNotQuantum is None:\n"
"            thermalNotQuantum = ( self.potentialScaler.currentTemperature > 0.0 )\n"
"        falseVacuumArray = FieldDictionaryToArray( falseVacuum[ \"FieldValues\" ] )\n"
"        trueVacuumArray = FieldDictionaryToArray( trueVacuum[ \"FieldValues\" ] )\n"
"        arrayOfArrays = numpy.array( [ trueVacuumArray, falseVacuumArray ] )\n"
"        configurationDifferenceLength = math.sqrt( numpy.sum( ( trueVacuumArray\n"
"                                                                - falseVacuumArray )**2 ) )\n"
"        tunnelPathPoints = int( configurationDifferenceLength / self.tunnelPathResolutionInGev )\n"
"        print( \"tunnelingCalculator.npoints being given \" + str( tunnelPathPoints ) )\n"
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
"            \n"
"        tunnelingCalculator = CTPD.fullTunneling( phi = arrayOfArrays,\n"
"                                                  V = self.potentialScaler.PotentialFromMatrix,\n"
"                                                  dV = self.potentialScaler.GradientFromMatrix,\n"
"                                                  alpha = numberOfSpaceTimeDimensionsForTunneling,\n"
"                                                  npoints = tunnelPathPoints,\n"
"                                                  quickTunneling = False )\n"
"\n"
"        if deformPath:\n"
"            tunnelingCalculator.run( **keywordArguments )\n"
"        else:\n"
"            tunnelingCalculator.tunnel1D()\n"
"\n"
"        calculatedAction = tunnelingCalculator.findAction()\n"
"\n"
"        if thermalNotQuantum:\n"
"            self.currentThermalAction = calculatedAction\n"
"        else:\n"
"            self.currentQuantumAction = calculatedAction\n"
"\n"
"        return calculatedAction\n"
"\n"
"\n"
"    def SetUpDirectQuantumTunneling( self, \n"
"                                     falseVacuum = None,\n"
"                                     trueVacuum = None ):\n"
"        if ( falseVacuum is None ):\n"
"            falseVacuum = self.dsbVacuum\n"
"        if ( trueVacuum is None ):\n"
"            trueVacuum = self.panicVacuum\n"
"        falseVacuumArray = FieldDictionaryToArray( falseVacuum[ \"FieldValues\" ] )\n"
"        trueVacuumArray = FieldDictionaryToArray( trueVacuum[ \"FieldValues\" ] )\n"
"        arrayOfArrays = numpy.array( [ trueVacuumArray, falseVacuumArray ] )\n"
"        configurationDifferenceLength = math.sqrt( numpy.sum( ( trueVacuumArray\n"
"                                                                - falseVacuumArray )**2 ) )\n"
"        print( \"quantumTunneler.npoints being given \"\n"
"               + str( int( configurationDifferenceLength / self.tunnelPathResolutionInGev ) ) )\n"
"        print( \"falseVacuumArray = \" + str( falseVacuumArray ) )\n"
"        print( \"trueVacuumArray = \" + str( trueVacuumArray ) )\n"
"        self.quantumTunneler = CTPD.fullTunneling( phi = arrayOfArrays,\n"
"                                                   V = self.potentialScaler.PotentialFromMatrix,\n"
"                                                   dV = self.potentialScaler.GradientFromMatrix,\n"
"                                                   alpha = 3,\n"
"                                                   npoints = int( configurationDifferenceLength / self.tunnelPathResolutionInGev ),\n"
"                                                   quickTunneling = False )\n"
"        self.quantumTunneler.tunnel1D( xtol = 1e-4, phitol = 1e-6 )\n"
"        self.quantumTunneler.fixEnd = False\n"
"\n"
"\n"
"    def CalculateQuantumAction( self ):\n"
"        self.previousQuantumAction = self.currentQuantumAction\n"
"        self.quantumTunneler.tunnel1D()\n"
"        self.currentQuantumAction = self.quantumTunneler.findAction()\n"
"        return self.currentQuantumAction\n"
"\n"
"\n"
"    def CalculateThermalAction( self ):\n"
"        self.previousThermalAction = self.currentThermalAction\n"
"        self.thermalTunneler.tunnel1D()\n"
"        self.currentThermalAction = self.thermalTunneler.findAction()\n"
"        return self.currentThermalAction\n"
"\n"
"\n"
"    def AllowedRunningTimeExceeded( self ):\n"
"        print( \"AllowedRunningTimeExceeded called. Running time so far = \"\n"
"               + str( time.clock() - self.startingTime )\n"
"               + \" seconds. Total allowed = \" + str( self.allowedRunningTime )\n"
"               + \" seconds. Allowed end time = \" + str( self.startingTime + self.allowedRunningTime )\n"
"               + \" seconds. Current time = \" + str( time.clock() )\n"
"               + \" seconds. About to return \" + str( time.clock() > ( self.startingTime + self.allowedRunningTime ) )\n"
"               + \".\" )\n"
"        return ( time.clock() > ( self.startingTime + self.allowedRunningTime ) )\n"
"\n"
"\n"
"\n"
"    def QuantumTunnelingTimeInInverseGev( self, currentAction ):\n"
"        if ( currentAction < 0.0 ):\n"
"            print( \"Warning! Quantum tunneling time for a negative action (\"\n"
"                   + str( currentAction )\n"
"                   + \") asked for, which may indicate a problem (or may not,\"\n"
"                   + \" but it would be a very thin potential barrier).\"\n"
"                   + \" Returning a tunneling time for an action of 0.\" )\n"
"            currentAction = 0.0\n"
"        elif ( currentAction > exponentCutOff ):\n"
"            print( \"Warning! Quantum tunneling time for a very large action (\"\n"
"                   + str( currentAction )\n"
"                   + \") asked for, which would probably cause an overflow in\"\n"
"                   + \" the evaluation of the exponent.\"\n"
"                   + \" Returning a tunneling time for an action of \"\n"
"                   + str( exponentCutOff ) + \".\" )\n"
"            currentAction = exponentCutOff\n"
"        return ( math.exp( 0.25 * currentAction ) / self.fourthRootOfSolitonicFactorA )\n"
"\n"
"\n"
"    def SetTemperature( self, currentTemperature ):\n"
"        self.potentialScaler.temperatureValue = currentTemperature\n"
"        if ( currentTemperature <= 0.0 ):\n"
"            self.thermalActionThreshold = -exponentCutOff\n"
"        else:\n"
"# Based on correspondence with Alexander Kusenko:\n"
"# Threshold Gamma / volume ~ T^8 / (N M_Planck^4) where N is the number of\n"
"# horizons.\n"
"# He took N = (10^14)^3. Assuming that the solitonic factor A is roughly T^4,\n"
"# then A exp( -S_{3-d,threshold} ) = 10^-42 T^8 / M_Planck^4\n"
"# so S_{3-d,threshold} = -T * ( -42 ln( 10 ) + 4 ln( T / M_Planck ) )\n"
"# = T * ( 42 ln( 10 ) + 4 ln( M_Planck/ T ) ).\n"
"            self.thermalActionThreshold = ( currentTemperature\n"
"                                            * ( lnOfNumberOfHorizons\n"
"                                                + ( 4.0 * math.log( reducedPlanckMass / currentTemperature ) ) ) )\n"
"            print( \"SetTemperature( \" + str( currentTemperature ) + \" ) called.\" )\n"
"            print( \"lnOfNumberOfHorizons = \" + str( lnOfNumberOfHorizons ) )\n"
"            print( \"reducedPlanckMass = \" + str( reducedPlanckMass ) )\n"
"            print( \"( reducedPlanckMass / currentTemperature ) = \" + str( reducedPlanckMass / currentTemperature ) )\n"
"            print( \"math.log( reducedPlanckMass / currentTemperature ) = \" + str( math.log( reducedPlanckMass / currentTemperature ) ) )\n"
"            print( \"self.thermalActionThreshold = \" + str( self.thermalActionThreshold ) )\n"
"\n"
"\n"
"    def DeformQuantumTunnelingPath( self ):\n"
"# maxiter is being passed to a CTPD.Deformation.deform( **deformationParams ) call.\n"
"        self.quantumTunneler.deform( nb = ( 2 * numberOfFields ), kb = 4, maxiter = 10 )\n"
"        self.quantumTunneler.tunnel1D()\n"
"\n"
"\n"
"    def DeformThermalTunnelingPath( self ):\n"
"# maxiter is being passed to a CTPD.Deformation.deform( **deformationParams ) call.\n"
"        self.thermalTunneler.deform( nb = ( 2 * numberOfFields ), kb = 4, maxiter = 10 )\n"
"        self.thermalTunneler.tunnel1D()\n"
"\n"
"\n"
"    def SetUpDirectThermalTunneling( self, \n"
"                                     falseVacuum = None,\n"
"                                     trueVacuum = None ):\n"
"        if ( falseVacuum is None ):\n"
"            falseVacuum = self.dsbVacuum\n"
"        if ( trueVacuum is None ):\n"
"            trueVacuum = self.panicVacuum\n"
"        falseVacuumArray = FieldDictionaryToArray( falseVacuum[ \"FieldValues\" ] )\n"
"        trueVacuumArray = FieldDictionaryToArray( trueVacuum[ \"FieldValues\" ] )\n"
"        arrayOfArrays = numpy.array( [ trueVacuumArray, falseVacuumArray ] )\n"
"        configurationDifferenceLength = math.sqrt( numpy.sum( ( trueVacuumArray\n"
"                                                                - falseVacuumArray )**2 ) )\n"
"        print( \"thermalTunneler.npoints being given \"\n"
"               + str( int( configurationDifferenceLength / self.tunnelPathResolutionInGev ) ) )\n"
"        self.thermalTunneler = CTPD.fullTunneling( phi = arrayOfArrays,\n"
"                                                   V = self.potentialScaler.PotentialFromMatrix,\n"
"                                                   dV = self.potentialScaler.GradientFromMatrix,\n"
"                                                   alpha = 2,\n"
"                                                   npoints = int( configurationDifferenceLength / self.tunnelPathResolutionInGev ),\n"
"                                                   quickTunneling = False )\n"
"        self.thermalTunneler.tunnel1D( xtol = 1e-4, phitol = 1e-6 )\n"
"        self.thermalTunneler.fixEnd = False\n"
"\n"
"\n"
"    def CalculateThermalAction( self ):\n"
"        self.previousThermalAction = self.currentThermalAction\n"
"        self.thermalTunneler.tunnel1D()\n"
"        self.currentThermalAction = self.thermalTunneler.findAction()\n"
"        return self.currentThermalAction\n"
"\n"
"    def ThermalTunnelingSurvivalProbability( self, exclusionTemperature, currentAction ):\n"
"        print( \"ThermalTunnelingSurvivalProbability( exclusionTemperature = \"\n"
"               + str( exclusionTemperature )\n"
"               + \", currentAction = \"\n"
"               + str( currentAction )\n"
"               + \" ) called.\" )\n"
"        actionOverTemperature = ( currentAction / exclusionTemperature )\n"
"        if ( actionOverTemperature > exponentCutOff ):\n"
"            print( \"Warning! Thermal tunneling survival probability for a very\"\n"
"                   \" large ratio (\"\n"
"                   + str( actionOverTemperature )\n"
"                   + \") of 3-dimensional action (\"\n"
"                   + str( currentAction )\n"
"                   + \") to temperature (\"\n"
"                   + str( exclusionTemperature )\n"
"                   + \") asked for, which would probably cause an overflow in\"\n"
"                   + \" the evaluation of the exponent.\"\n"
"                   + \" Returning a survival probability of 0.\" )\n"
"            return 0.0\n"
"        else:\n"
"            decayWidth = ( ( numberOfHorizons\n"
"                             * reducedPlanckMass**3\n"
"                             * math.exp( -actionOverTemperature ) )\n"
"                           / ( exclusionTemperature**2 ) )\n"
"            timeAtTemperature = ( reducedPlanckMass / ( exclusionTemperature**2 ) )\n"
"            exponentiationFactor = ( -decayWidth * timeAtTemperature )\n"
"            print( \"decayWidth = \" + str( decayWidth ) )\n"
"            print( \"timeAtTemperature = \" + str( timeAtTemperature ) )\n"
"            print( \"exponentiationFactor = \" + str( exponentiationFactor ) )\n"
"            if ( exponentiationFactor < -exponentCutOff ):\n"
"                exponentiationFactor = -exponentCutOff\n"
"                print( \"Warning! Thermal decay width times time at temperature\"\n"
"                       + \" is so large and negative that an overflow error\"\n"
"                       + \" would occur when exponentiating to get the survival\"\n"
"                       + \" probability, so capping it at \"\n"
"                       + str( exponentiationFactor )\n"
"                       + \".\" )\n"
"            if ( exponentiationFactor > exponentCutOff ):\n"
"                exponentiationFactor = exponentCutOff\n"
"                print( \"Warning! Thermal decay width times time at temperature\"\n"
"                       + \" is so large and positive that an overflow error\"\n"
"                       + \" would occur when exponentiating to get the survival\"\n"
"                       + \" probability, so capping it at \"\n"
"                       + str( exponentiationFactor )\n"
"                       + \".\" )\n"
"            print( \"returning \" + str( math.exp( exponentiationFactor ) ) )\n"
"            return math.exp( exponentiationFactor )\n"
"\n"
"# End of Vevacious class\n"
"\n"














<< vevaciousVersionName << " = \"" << vevaciousVersionString << "\"\n"
<< resultsFilenameVariableName << " = \"" << resultsFilename << "\"\n"
"\n"
<< PotentialMinimizer::energyScale << " = " << sarahInterpreter.getSlhaScale()
<<                                                                         "\n"
<< saddleSplitNudges << " = [ ";
    for( std::vector< double >::iterator
         whichNudge( saddleNudgeList.begin() );
         saddleNudgeList.end() > whichNudge;
         ++whichNudge )
    {
      if( saddleNudgeList.begin() != whichNudge )
      {
        outputFile << ", ";
      }
      outputFile
      << "( " << *whichNudge << " / " << PotentialMinimizer::energyScale
      <<                                                                  " )";
    }
    outputFile << " ]\n"
"\n"
<< rollingToleranceVariableName << " = " << rollingTolerance << "\n"
<< pathToCosmotransitionsVariableName << " = \"" << pathToCosmotransitions
<<                                                                       "\"\n"
<< directLifetimeBoundVariableName << " = " << directLifetimeBound << "\n"
<< deformedLifetimeBoundVariableName << " = " << lifetimeThreshold << "\n";
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
          else if( 0 == stabilityVerdict.compare( "long-lived" ) )
          {
            stabiltyResult = 0.0;
          }
          else if( 0 == stabilityVerdict.compare( "short-lived" ) )
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
      << "BLOCK VEVACIOUSRESULTS # results from Vevacious version "
      << vevaciousVersionString << ", documented in " << vevaciousDocumentation
      << "\n"
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
"# (B.O'L. would like to apologize about the state of this file: ideally it\n"
"# would be a lot neater, and respect modern programming practices. However,\n"
"# it ended up easier to try to keep all this in one file. Maybe a later\n"
"# version will tidy it up properly, maybe splitting it into separate \n"
"# files...)\n"
"#\n"
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
"\n"
"effectivePotentialFunction = VPD."
<<                           PotentialMinimizer::loopCorrectedPotential << "\n"
"# one could replace VPD." << PotentialMinimizer::loopCorrectedPotential
<<                                                                         "\n"
"# with VPD." << PotentialMinimizer::treeLevelPotential << " for a\n"
"# tree-level analysis, for example.\n"
"\n"
"potentialAtVevOrigin = ( VPD." << PotentialMinimizer::energyScaleFourth
<<                                                                         "\n"
"                         * VPD." << PotentialMinimizer::functionFromDictionary
<<                                            "( effectivePotentialFunction,\n"
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
"minuitObject = minuit.Minuit( effectivePotentialFunction )\n"
"# There is a possibility of MINUIT trying to roll off to infinity but\n"
"# generating an overflow error in calculating the potential before throwing\n"
"# an exception. Here limits are set on the VEVs of one hundred times\n"
"# the energy scale, which is probably way beyond where the potential should\n"
"# be trusted anyway.\n"
"for vevVariable in minuitObject.values:\n"
"    minuitObject.limits[ vevVariable ] = ( -100.0, 100.0 )\n"
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
"                                                effectivePotentialFunction,\n"
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
"\n"
"# the global minimum of the tree-level potential is also recorded, in case\n"
"# loop corrections move around basins of attraction too far.\n"
"globalTreeMinimum = pointsToTry[ 0 ]\n"
"for vevKey in globalTreeMinimum.keys():\n"
"    globalTreeMinimum[ vevKey ] = 0.0\n"
"globalTreeMinimumDepth = VPD." << PotentialMinimizer::functionFromDictionary
<<                  "( VPD." << PotentialMinimizer::treeLevelPotential << ",\n"
"                                                        globalTreeMinimum )\n"
"for vevValueSet in " << pointsToTry << ":\n"
"    if ( VevsHaveCorrectSigns( vevValueSet ) ):\n"
"        TryToMinimize( vevValueSet )\n"
"        treeLevelDepth = VPD." << PotentialMinimizer::functionFromDictionary
<<                  "( VPD." << PotentialMinimizer::treeLevelPotential << ",\n"
"                                             vevValueSet )\n"
"        if ( treeLevelDepth < globalTreeMinimumDepth ):\n"
"            globalTreeMinimum = vevValueSet.copy()\n"
"            globalTreeMinimumDepth = treeLevelDepth\n"
"\n"
"for saddleSplitNudge in VPD." << saddleSplitNudges << ":\n"
"    if ( 0 < len( foundSaddles ) ):\n"
"        print( \"PyMinuit had to be nudged off \"\n"
"               + str( len( foundSaddles ) ) + \" saddle point(s).\" )\n"
"        " << pointsToTry << " = []\n"
"        for saddlePoint in foundSaddles:\n"
"            " << pointsToTry << ".append( DisplacePoint(\n"
"                                         saddlePoint[ 0 ][ \'vevValues\' ],\n"
"                                                       saddlePoint[ 1 ], \n"
"                                                       saddleSplitNudge ) )\n"
"            " << pointsToTry << ".append( DisplacePoint(\n"
"                                         saddlePoint[ 0 ][ \'vevValues\' ],\n"
"                                                       saddlePoint[ 1 ], \n"
"                                                      -saddleSplitNudge ) )\n"
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
"# The result is assumed long-lived metastable unless found otherwise.\n"
"stabilityVerdict = \"long-lived\"\n"
"givenInputAsArray = VPD." << PotentialMinimizer::vevDictionaryToArray
<<                     "( VPD." << PotentialMinimizer::inputVevsPoint << " )\n"
"minuitObject.values = VPD." << PotentialMinimizer::inputVevsPoint
<<                                                                  ".copy()\n"
"foundSaddles = []\n"
"try:\n"
"    minuitObject.migrad()\n"
"    candidateDescent = SteepestDescent( minuitObject.values )\n"
"    if ( 0.0 > candidateDescent[ 0 ] ):\n"
"        warningMessage = (\n"
"                     \"Input VEVs seem to correspond to a saddle point!\" )\n"
"        warningMessages.append( warningMessage )\n"
"        print( warningMessage )\n"
"except minuit.MinuitError as minuitError:\n"
"    warningMessage = ( \"PyMinuit had problems starting at input VEVs! \"\n"
"                     + VPD."
<<                        PotentialMinimizer::userVevsAsMathematica << "( VPD."
<<                                 PotentialMinimizer::inputVevsPoint << " )\n"
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
"# this is a check that there wasn't a weird expansion of the basin of\n"
"# attraction of the loop input minimum to include the tree minimum that\n"
"# should have rolled to an undesired loop minimum that moved, with its\n"
"# basin of attraction, quite far away.\n"
"if ( ( rollingToleranceSquared * rolledInputLengthSquared )\n"
"     > distanceSquaredFromInputToGlobalMinimum ):\n"
"    globalTreeMinimumAsArray = VPD."
<<        PotentialMinimizer::vevDictionaryToArray << "( globalTreeMinimum )\n"
"    if ( ( rollingToleranceSquared * rolledInputLengthSquared )\n"
"         < numpy.sum( ( globalTreeMinimumAsArray\n"
"                        - rolledInputAsArray )**2 ) ):\n"
"        warningMessage = ( \"Initially judged stable at 1-loop level, but\"\n"
"                          + \" (possibly) metastable at tree level!\"\n"
"                          + \" Adding scaled tree-level global minimum as\"\n"
"                          + \" an extra PyMinuit starting point.\" )\n"
"        warningMessages.append( warningMessage )\n"
"        print( warningMessage )\n"
"        for vevKey in globalTreeMinimum.keys():\n"
"            globalTreeMinimum[ vevKey ] = ( 2.0\n"
"                                            * globalTreeMinimum[ vevKey ] )\n"
"        minuitObject.values = globalTreeMinimum.copy()\n"
"        try:\n"
"            minuitObject.migrad()\n"
"        except minuit.MinuitError as minuitError:\n"
"            warningMessage = ( \"PyMinuit had problems starting at\"\n"
"                               + \" doubled VEVs of tree-level global\"\n"
"                               + \" minimum! \"\n"
"                               + VPD."
<<                        PotentialMinimizer::userVevsAsMathematica << "( VPD."
<<                                 PotentialMinimizer::inputVevsPoint << " )\n"
"                               + \" [minuit.MinuitError: \"\n"
"                               + str( minuitError )\n"
"                               + \"]. PyMinuit stopped at \"\n"
"                          + VPD." << PotentialMinimizer::userVevsAsMathematica
<<                                                  "( minuitObject.values )\n"
"                               + \" with relative depth \"\n"
"               + str( ( VPD." << PotentialMinimizer::energyScaleFourth << "\n"
"                        * VPD."
<<                       PotentialMinimizer::functionFromDictionary << "( VPD."
<<                          PotentialMinimizer::loopCorrectedPotential << ",\n"
"                                                    minuitObject.values ) )\n"
"                      - potentialAtVevOrigin )\n"
"                               + \" at one-loop level and \"\n"
"                 + str( VPD." << PotentialMinimizer::energyScaleFourth << "\n"
"                        * VPD."
<<                       PotentialMinimizer::functionFromDictionary << "( VPD."
<<                              PotentialMinimizer::treeLevelPotential << ",\n"
"                                                    minuitObject.values ) )\n"
"                               + \" at tree level.\"\n"
"               + \" Minuit's estimate of how much deeper it should go is \"\n"
"                 + str( VPD." << PotentialMinimizer::energyScaleFourth << "\n"
"                        * minuitObject.edm ) )\n"
"            warningMessages.append( warningMessage )\n"
"            print( warningMessage )\n"
"        globalMinimumPointAsDictionary = minuitObject.values.copy()\n"
"        globalMinimumPointAsArray = VPD."
<<      PotentialMinimizer::vevDictionaryToArray << "( minuitObject.values )\n"
"        globalMinimumDepthValue = minuitObject.fval\n"
"        distanceSquaredFromInputToGlobalMinimum = numpy.sum( (\n"
"                      globalMinimumPointAsArray - rolledInputAsArray )**2 )\n"
"        actionNeedsToBeCalculated = True\n"
"        print( \"scaled point rolled to \"\n"
"                       + VPD." << PotentialMinimizer::userVevsAsMathematica
<<                                                  "( minuitObject.values )\n"
"               + \" [distanceSquaredFromInputToGlobalMinimum = \"\n"
"               + str( distanceSquaredFromInputToGlobalMinimum )\n"
"         + \" ; ( rollingToleranceSquared * rolledInputLengthSquared ) = \"\n"
"            + str( ( rollingToleranceSquared * rolledInputLengthSquared ) )\n"
"                       + \"] with relative depth \"\n"
"               + str( ( VPD." << PotentialMinimizer::energyScaleFourth << "\n"
"                        * VPD."
<<                       PotentialMinimizer::functionFromDictionary << "( VPD."
<<                          PotentialMinimizer::loopCorrectedPotential << ",\n"
"                                                    minuitObject.values ) )\n"
"                      - potentialAtVevOrigin ) )\n"
"        if ( globalMinimumDepthValue >= rolledInputDepth ):\n"
"            globalMinimumPointAsDictionary = rolledInputAsDictionary\n"
"            globalMinimumPointAsArray = rolledInputAsArray\n"
"            globalMinimumDepthValue = rolledInputDepth\n"
"            distanceSquaredFromInputToGlobalMinimum = 0.0\n"
"            actionNeedsToBeCalculated = False\n"
"\n"
"\n"
"# The decay width per unit volume, following Sidney Coleman, is of the form\n"
"# A exp( -[bounce action] ), where A is a solitonic solution that is\n"
"# usually actually just estimated on dimensional grounds, as since the\n"
"# action needs to be about 400 for a tunneling time of about 14 billion\n"
"# years, so getting the action right to a per-cent is more important than\n"
"# getting the A factor right within a factor of 3! Hence we just take A to\n"
"# be the scale from the SLHA file to the fourth power, which should be\n"
"# about the same size as the largest dimensionful Lagrangian parameter\n"
"# entering the effective potential (to stabilize loop corrections), so\n"
"# coincidentally should be the same order of magnitude as one would expect\n"
"# the solitonic solutions to be. The user can modify this here, if they\n"
"# feel that it should be something else, e.g. 0.1 times this value, or that\n"
"# instead it should be of the order of the electroweak VEV of 246 GeV.\n"
"# fourthRootOfSolitonicFactorA should be in units of GeV (as it is\n"
"# A^(1/4)).\n"
"fourthRootOfSolitonicFactorA = VPD." << PotentialMinimizer::energyScaleFourth
<<                                                                         "\n"
"\n"
"# The age of the known Universe is pretty close to exactly (10^(41))/GeV.\n"
"ageOfTheKnownUniverseInInverseGev = 1.0E+41\n"
"\n"
"# The tunneling time calculation action thresholds are initially turned\n"
"# off, but if positive lifetime thresholds were given, the action\n"
"# thresholds get set correctly.\n"
"directActionThreshold = -1\n"
"deformedActionThreshold = -1\n"
"# The tunneling time should be something like\n"
"# (decay width / unit volume)^(1/4), which explains the factor of 4.0 in\n"
"# the following thresholds on the bounce action.\n"
"# (This is from:\n"
"# Gamma_threshold / unit volume = t_threshold^(-4) = A exp( -B_threshold )\n"
"# -B_threshold = ln( 1 / ( A t_threshold^4 ) )\n"
"# B_threshold = 4 ln( A^(1/4) t_threshold )\n"
"# hence the factor of 4 for the threshold actions.)\n"
"if ( 0.0 < VPD." << directLifetimeBoundVariableName << " ):\n"
"    directActionThreshold = ( 4.0 * math.log( VPD."
<<                                      directLifetimeBoundVariableName << "\n"
"                                        * ageOfTheKnownUniverseInInverseGev\n"
"                                         * fourthRootOfSolitonicFactorA ) )\n"
"if ( 0.0 < VPD." << deformedLifetimeBoundVariableName << " ):\n"
"    deformedActionThreshold = ( 4.0 * math.log( VPD."
<<                                    deformedLifetimeBoundVariableName << "\n"
"                                        * ageOfTheKnownUniverseInInverseGev\n"
"                                         * fourthRootOfSolitonicFactorA ) )\n"
"\n"
"\n"
"# If the input vacuum is the global minimum, actionValue is set to -1.0.\n"
"# Non-negative values of actionValue indicate the current upper bound on\n"
"# the action after the last approximation. It starts stupidly high\n"
"# (the current age of the Universe corresponds to an action of about 400,\n"
"# for an A factor of (100 GeV)^4) so that the first requested bounding\n"
"# estimate will be calculated.\n"
"if ( ( rollingToleranceSquared * rolledInputLengthSquared )\n"
"     >= distanceSquaredFromInputToGlobalMinimum ):\n"
"    stabilityVerdict = \"stable\"\n"
"    actionValue = -1.0\n"
"    tunnelingTime = -1.0\n"
"    actionType = \"unnecessary\"\n"
"    actionNeedsToBeCalculated = False\n"
"\n"
"# The resolution of the tunneling path needs to be set\n"
"# (low-ish by default for speed):\n"
"tunnelingResolution = 20\n"
"\n"
"if ( actionNeedsToBeCalculated\n"
"     and ( ( 0.0 < directActionThreshold )\n"
"           or ( 0.0 < deformedActionThreshold ) ) ):\n"
"    firstStepPoint = ( ( rolledInputAsArray\n"
"                         * ( 1.0 - ( 1.0 / tunnelingResolution ) ) )\n"
"                       + ( globalMinimumPointAsArray\n"
"                           * ( 1.0 / tunnelingResolution ) ) )\n"
"    firstStepDepth = VPD." << PotentialMinimizer::functionFromArray
<<                                            "( effectivePotentialFunction,\n"
"                                                firstStepPoint )\n"
"    if ( rolledInputDepth >= firstStepDepth ):\n"
"        actionType = \"barrier_smaller_than_resolution\"\n"
"        actionValue = 0.0\n"
"        tunnelingTime = 0.0\n"
"        stabilityVerdict = \"short-lived\"\n"
"        actionNeedsToBeCalculated = False\n"
"        warningMessage = ( \"Energy barrier from input VEVs to global\"\n"
"                        + \" minimum thinner than resolution of tunneling\"\n"
"                           + \" path!\" )\n"
"        warningMessages.append( warningMessage )\n"
"        print( warningMessage )\n"
"\n"
"if ( actionNeedsToBeCalculated\n"
"     and ( ( 0.0 < directActionThreshold )\n"
"           or ( 0.0 < deformedActionThreshold ) ) ):\n"
"    import sys\n"
"    sys.path.append( VPD." << pathToCosmotransitionsVariableName << " )\n"
"    import pathDeformation as CPD\n"
"\n"
"    arrayOfArrays = numpy.array( [ globalMinimumPointAsArray.copy(),\n"
"                                   rolledInputAsArray.copy() ] )\n"
"\n"
"    def PotentialFromArray( pointAsArray ):\n"
"        return VPD." << PotentialMinimizer::functionFromArray
<<                                            "( effectivePotentialFunction,\n"
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
"    if ( ( 0.0 < directActionThreshold )\n"
"         and ( actionValue > directActionThreshold ) ):\n"
"        quickTunneler = CPD.fullTunneling( phi = arrayOfArrays,\n"
"                                           V = PotentialFromMatrix,\n"
"                                           dV = GradientFromMatrix,\n"
"                                           alpha = 3,\n"
"                                           npoints = tunnelingResolution,\n"
"                                           quickTunneling = False )\n"
"        quickTunneler.tunnel1D( xtol = 1e-4, phitol = 1e-6 )\n"
"        actionValue = quickTunneler.findAction()\n"
"        actionType = \"direct_path_bound\"\n"
"        if( actionValue < directActionThreshold ):\n"
"            stabilityVerdict = \"short-lived\"\n"
"            actionNeedsToBeCalculated = False\n"
"\n"
"    if ( actionNeedsToBeCalculated\n"
"         and ( 0.0 < deformedActionThreshold )\n"
"         and ( actionValue > deformedActionThreshold ) ):\n"
"        fullTunneler = CPD.fullTunneling( phi = arrayOfArrays,\n"
"                                           V = PotentialFromMatrix,\n"
"                                           dV = GradientFromMatrix,\n"
"                                           alpha = 3,\n"
"                                           npoints = tunnelingResolution,\n"
"                                           quickTunneling = False )\n"
"# setting a maximum of 20 path deformation iterations before giving up on\n"
"# finding the optimal path may seem defeatist, but I have rarely seen it\n"
"# converge if it hasn't within the first few iterations. an action is still\n"
"# calculated, though it may not be the minimum action possible.\n"
"        fullTunneler.run( maxiter = 20 )\n"
"        actionValue = fullTunneler.findAction()\n"
"        actionType = \"full_deformed_path\"\n"
"        if ( actionValue < deformedActionThreshold ):\n"
"            stabilityVerdict = \"short-lived\"\n"
"            actionNeedsToBeCalculated = False\n"
"\n"
"# No matter if there were serious errors or not, an output file is written:\n"
"outputFile = open( VPD." << resultsFilenameVariableName << ", \"w\" )\n"
"outputFile.write( \"<Vevacious_result>\\n\"\n"
"                  + \"  <reference version=\\\"" << vevaciousVersionString
<<           "\\\" citation=\\\"" << vevaciousDocumentation << "\\\" />\\n\"\n"
"             + \"  <stability> \" + stabilityVerdict + \" </stability>\\n\"\n"
"                  + \"  <global_minimum   relative_depth=\\\"\"\n"
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
"# The tunneling time should be something like\n"
"# (decay width / unit volume)^(1/4), which explains the factor of 0.25 in\n"
"# the exponential here.\n"
"    tunnelingTime = ( math.exp( 0.25 * actionValue )\n"
"                      / ( ageOfTheKnownUniverseInInverseGev\n"
"                          * fourthRootOfSolitonicFactorA ) )\n"
"    if ( 1000000.0 < tunnelingTime ):\n"
"        tunnelingTime = 1000000.0\n"
"# The tunneling time is given in units of the current age of the known\n"
"# Universe, and is capped at one million.\n"
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
