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
  std::string const VevaciousRunner::vevaciousVersionString( "1.0.7" );
  std::string const
  VevaciousRunner::vevaciousVersionName( "vevaciousVersion" );
  std::string const
  VevaciousRunner::vevaciousDocumentation( "arXiv:1307.1477 (hep-ph)" );
  std::string const
  VevaciousRunner::pointsToTry( "pointsToTry" );
  std::string const
  VevaciousRunner::resultsFilenameVariableName( "outputFile" );
  std::string const
  VevaciousRunner::saddleSplitNudges( "saddleSplitNudges" );
  std::string const
  VevaciousRunner::rollingToleranceVariableName( "rollingTolerance" );
  std::string const
  VevaciousRunner::quarticActionBoundVariableName( "quarticActionBound" );
  std::string const VevaciousRunner::pathToCosmotransitionsVariableName(
                                                    "pathToCosmotransitions" );
  std::string const
  VevaciousRunner::directLifetimeBoundVariableName( "directLifetimeBound" );
  std::string const VevaciousRunner::deformedLifetimeBoundVariableName(
                                                     "deformedLifetimeBound" );
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
    directLifetimeBound( 0.1 ),
    deformedLifetimeBound( directLifetimeBound )
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
    directLifetimeBound( BOL::StringParser::stringToDouble(
                                         argumentParser.fromTag( "direct_time",
                                                                 "0.1" ) ) ),
    deformedLifetimeBound( BOL::StringParser::stringToDouble(
                                       argumentParser.fromTag( "deformed_time",
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
<< deformedLifetimeBoundVariableName << " = " << deformedLifetimeBound << "\n"
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
"globalTreeMinimumDepth = VPD." << PotentialMinimizer::functionFromDictionary
<<                  "( VPD." << PotentialMinimizer::treeLevelPotential << ",\n"
"\n                                                    globalTreeMinimum )\n"
"for vevValueSet in " << pointsToTry << ":\n"
"    if ( VevsHaveCorrectSigns( vevValueSet ) ):\n"
"        TryToMinimize( vevValueSet )\n"
"treeLevelDepth = VPD." << PotentialMinimizer::functionFromDictionary
<<                  "( VPD." << PotentialMinimizer::treeLevelPotential << ",\n"
"                                             vevValueSet )\n"
"if ( treeLevelDepth < globalTreeMinimum ):\n"
"    globalTreeMinimum = vevValueSet.copy()\n"
"    globalTreeMinimumDepth = treeLevelDepth\n"
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
"# fourthRootOfSolitonicFactorA should be in units of GeV (as it is A^(1/4)\n"
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
"        quickTunneler = CPD.fullTunneling( V = PotentialFromMatrix,\n"
"                                           dV = GradientFromMatrix,\n"
"                                           phi = arrayOfArrays,\n"
"                                           quickTunneling = False,\n"
"                                           npoints = tunnelingResolution )\n"
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
"        fullTunneler = CPD.fullTunneling( V = PotentialFromMatrix,\n"
"                                          dV = GradientFromMatrix,\n"
"                                          phi = arrayOfArrays,\n"
"                                          quickTunneling = False,\n"
"                                          npoints = tunnelingResolution )\n"
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
