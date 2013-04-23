/*
 * PotentialMinimizer.cpp
 *
 *  Created on: Oct 5, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of Vevacious, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.Vevacious.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "../include/PotentialMinimizer.hpp"

namespace Vevacious
{
  std::string const PotentialMinimizer::namesOfVevs( "namesOfVevs" );
  std::string const PotentialMinimizer::internalVevNamesToUserNames(
                                               "internalVevNamesToUserNames" );
  std::string const
  PotentialMinimizer::vevsTakenPositive( "vevsTakenPositive" );
  std::string const PotentialMinimizer::energyScale( "energyScale" );
  std::string const
  PotentialMinimizer::energyScaleFourth( "energyScaleFourth" );
  std::string const
  PotentialMinimizer::userVevsAsMathematica( "UserVevsAsMathematica" );
  std::string const
  PotentialMinimizer::userVevsAsXml( "UserVevsAsXml" );
  std::string const
  PotentialMinimizer::vevDictionaryToArray( "VevDictionaryToArray" );
  std::string const PotentialMinimizer::vevOrigin( "vevOrigin" );
  std::string const PotentialMinimizer::inputVevsPoint( "inputVevsPoint" );
  std::string const
  PotentialMinimizer::functionFromDictionary( "FunctionFromDictionary" );
  std::string const
  PotentialMinimizer::functionFromArray( "FunctionFromArray" );
  std::string const
  PotentialMinimizer::loopCorrectedPotential( "LoopCorrectedPotential" );
  std::string const
  PotentialMinimizer::treeLevelPotential( "TreeLevelPotential" );

  std::string const
  PotentialMinimizer::inverseScaleSquared( "inverseScaleSquared" );
  std::string const
  PotentialMinimizer::inverseScaleFourthed( "inverseScaleFourthed" );
  std::string const PotentialMinimizer::loopFactor( "loopFactor" );
  std::string const PotentialMinimizer::overallFactorAttributeName( "factor" );
  std::string const
  PotentialMinimizer::subtractionConstantAttributeName( "constant" );
  std::string const
  PotentialMinimizer::defaultSubtractionConstantString( "1.5" );
  std::string const
  PotentialMinimizer::massSquaredEigenvalueNameBase( "MassesSquared" );
  std::string const
  PotentialMinimizer::massCorrectionNameBase( "MassCorrection" );
  std::string const PotentialMinimizer::polynomialPartFunctionName(
                                                "PolynomimalPartOfPotential" );
  std::string const
  PotentialMinimizer::massSquaredCorrections( "MassSquaredCorrections" );
  std::string const
  PotentialMinimizer::loopCorrectionFunctionName( "LoopCorrectionSum" );
  std::string const PotentialMinimizer::loopCorrectedPotentialFunctionName(
                                                    "LoopCorrectedPotential" );


  PotentialMinimizer::PotentialMinimizer() :
    sarahParser( false ),
    attributeFinder(),
    readSubtractionConstant( "" ),
    pythonCode(),
    polynomialPartOfPotential( "" ),
    factoredMassSquaredMatrices(),
    correctionFunctions(),
    stringParser()
  {
    pythonCode.precision( 17 );
    stringParser.precision( 17 );
    // HOM4PS2 operates at 17 digits of precision.
    // SARAH should prepare all the mass-squared matrices (even for fermions)
    // as XML elements. there is no set number of mass matrices.
  }

  PotentialMinimizer::~PotentialMinimizer()
  {
    sarahParser.closeFile();
  }


  void PotentialMinimizer::addMassSquaredMatrix(
                 std::map< std::string, std::string > const& factorAndConstant,
                                         std::string const& massSquaredMatrix )
  {
    attributeFinder
    = factorAndConstant.find( subtractionConstantAttributeName );
    if( factorAndConstant.end() != attributeFinder )
    {
      readSubtractionConstant.assign( BOL::StringParser::trimFromFrontAndBack(
                                                       attributeFinder->second,
                                                                    "\"\'" ) );
    }
    else
    {
      readSubtractionConstant.assign( defaultSubtractionConstantString );
    }
    attributeFinder = factorAndConstant.find( overallFactorAttributeName );
    if( factorAndConstant.end() != attributeFinder )
    {
      factoredMassSquaredMatrices.push_back( StringTriple() );
      factoredMassSquaredMatrices.back().first.first.assign(
              BOL::StringParser::trimFromFrontAndBack( attributeFinder->second,
                                                       "\"\'" ) );
      factoredMassSquaredMatrices.back().first.second.assign(
          readSubtractionConstant );
      factoredMassSquaredMatrices.back().second.assign( massSquaredMatrix );
      SarahInterpreter::removeNewlinesFrom(
                                   factoredMassSquaredMatrices.back().second );
    }
  }

  void PotentialMinimizer::prepareBasePython(
                                           SarahInterpreter& sarahInterpreter )
  {
    std::vector< char > const&
    vevsOrderedBySolution( sarahInterpreter.getVevsOrderedBySolutions() );
    pythonCode.clear();
    pythonCode.str( "" );
    pythonCode <<
"# Unfortunately PyMinuit is stupidly restricted to single-character\n"
"# variable names, so Vevacious internally renames the VEVs:\n"
"# " << sarahInterpreter.getHumanReadableVevNameMap() << "\n"
<< namesOfVevs << " = [ "
<<             sarahInterpreter.getInternalVevNamesAsQuotedCharList() << " ]\n"
<< internalVevNamesToUserNames << " = { ";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        pythonCode << ", ";
      }
      pythonCode << "\'" << *whichName << "\': \""
      << sarahInterpreter.getUserVevName( *whichName ) << "\"";
    }
    pythonCode << " }\n"
<< vevsTakenPositive << " = " << sarahInterpreter.getPositiveInternalVevs()
<<                                                                         "\n"
"\n"
<< energyScale << " = " << sarahInterpreter.getSlhaScale() << "\n"
<< energyScaleFourth << " = ( " << energyScale << " * " << energyScale << " * "
<<                                energyScale << " * " << energyScale << " )\n"
<< inverseScaleSquared << " = ( 1.0 / ( " << energyScale << " * "
<<                                                      energyScale << " ) )\n"
<< inverseScaleFourthed << " = ( 1.0 / " << energyScaleFourth << " )\n"
"\n"
"def " << userVevsAsMathematica << "( vevDictionary ):\n"
"    return \"{ ";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        pythonCode << ", ";
      }
      pythonCode << sarahInterpreter.getUserVevName( *whichName )
      << " -> ( \" + str( " << energyScale << " * vevDictionary[ \'"
      << *whichName << "\' ] ) + \" )";
    }
    pythonCode << " }\"\n"
"\n"
"def " << userVevsAsXml << "( vevDictionary ):\n"
"    return \"";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        pythonCode << " ";
      }
      pythonCode << sarahInterpreter.getUserVevName( *whichName )
      << "=\\\"\" + str( " << energyScale << " * vevDictionary[ \'"
      << *whichName << "\' ] ) + \"\\\"";
    }
    pythonCode << "\"\n"
"\n"
<< vevOrigin << " = { ";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        pythonCode << ", ";
      }
      pythonCode << "\'" << *whichName << "\': 0.0";
    }
    pythonCode << " }\n"
"\n"
<< inputVevsPoint << " = { ";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        pythonCode << ", ";
      }
      pythonCode << "\'" << *whichName << "\': ( "
      << sarahInterpreter.getInputVevValue( *whichName ) << " / "
      << energyScale << " )";
    }
    pythonCode << " }\n"
"\n"
"def " << vevDictionaryToArray << "( vevDictionary ):\n"
"    return numpy.array( [ ";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        pythonCode << ", ";
      }
      pythonCode << "vevDictionary[ \'" << *whichName << "\' ]";
    }
    pythonCode << " ] )\n"
"\n"
"def " << functionFromDictionary
<<                              "( FunctionFromArguments, vevDictionary ):\n"
"    return FunctionFromArguments( \n";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      if( vevsOrderedBySolution.begin() != whichName )
      {
        pythonCode << ", ";
      }
      pythonCode
      << *whichName << " = vevDictionary[ \'" << *whichName << "\' ]";
    }
    pythonCode << " )\n"
"\n"
"def " << functionFromArray << "( FunctionFromArguments, vevArray ):\n"
"    return FunctionFromArguments( \n";
    for( unsigned int whichVev( 0 );
         vevsOrderedBySolution.size() > whichVev;
         ++whichVev )
    {
      if( 0 < whichVev )
      {
        pythonCode << ", ";
      }
      pythonCode << vevsOrderedBySolution[ whichVev ] << " = vevArray[ "
      << whichVev << " ]";
    }
    pythonCode << " )\n"
"\n"
<< loopFactor << " = ( 1.0 / ( 64.0 * math.pi * math.pi ) )\n"
"\n"
"# The loop corrections are scaled to the energy scale, since the masses are\n"
"# scaled to the energy scale:\n"
"def " << massSquaredCorrections << "( massSquaredValues,\n"
"                                      overallFactor,\n"
"                                      subtractionConstant ):\n"
"    summedCorrection = 0.0\n"
"    for massSquaredValue in massSquaredValues:\n"
"        massSquaredMagnitude = abs( massSquaredValue )\n"
"        if( ( 1.0E-6 * " << inverseScaleSquared
<<                                             " ) < massSquaredMagnitude ):\n"
"            summedCorrection += ( massSquaredMagnitude\n"
"                                  * massSquaredMagnitude \n"
"                                  * ( math.log( massSquaredMagnitude )\n"
"                                      - subtractionConstant ) )\n"
"    return ( overallFactor * summedCorrection )\n"
"\n";
  }

  void PotentialMinimizer::prepareLoopCorrections(
                                           SarahInterpreter& sarahInterpreter )
  {
    correctionFunctions.clear();
    std::string massSquaredMatrix( "" );
    std::string overallFactor( "" );
    std::string subtractionConstant( "" );
    std::string massSquaredFunctionName( "" );
    BOL::VectorlikeArray< std::string > massTerms;
    for( std::vector< StringTriple >::iterator
         whichTriple( factoredMassSquaredMatrices.begin() );
         factoredMassSquaredMatrices.end() > whichTriple;
         ++whichTriple )
    {
      massTerms.clearEntries();
      BOL::StringParser::parseByChar( sarahInterpreter( whichTriple->second ),
                                      massTerms,
                                      ',' );
      stringParser.clear();
      stringParser.str( "" );
      stringParser
      << massSquaredEigenvalueNameBase << ( correctionFunctions.size() + 1 );
      massSquaredFunctionName.assign( stringParser.str() );
      pythonCode <<
"# The masses-squared from this function are scaled to the energy scale:\n"
"def " << massSquaredFunctionName << argumentsString << ":\n"
"    massSquaredElementArray = numpy.array( [ ";
      for( int whichTerm( 0 );
           massTerms.getSize() > whichTerm;
           ++whichTerm )
      {
        if( 0 < whichTerm )
        {
          pythonCode << ", ";
        }
        pythonCode << " ( "
        << BOL::StringParser::trimFromFrontAndBack( massTerms[ whichTerm ],
                                 BOL::StringParser::whitespaceAndNewlineChars )
        << " )";
      }
      pythonCode << " ] )\n"
"    massSquaredElementArray *= " << inverseScaleSquared << "\n"
"    massSquaredMatrix = massSquaredElementArray.reshape( math.sqrt(\n"
"                                       massSquaredElementArray.size ), -1 )\n"
"    return numpy.linalg.eigvalsh( massSquaredMatrix )\n"
"\n";
      stringParser.clear();
      stringParser.str( "" );
      stringParser
      << massCorrectionNameBase << ( correctionFunctions.size() + 1 );
      // the corrections are named "MassCorrection1(" + arguments + ")",
      // "MassCorrection2(" + arguments + ")" & so on.
      correctionFunctions.push_back( stringParser.str() );
      // all the corrections have to be noted for later.
      pythonCode <<
"def " << stringParser.str() << argumentsString << ":\n"
"    return " << massSquaredCorrections << "( " << massSquaredFunctionName
<< argumentsString << ", ( " << whichTriple->first.first << " ) , ( "
<< whichTriple->first.second << " ) )\n"
"\n";
    }
    pythonCode <<
"def " << loopCorrectionFunctionName << argumentsString << ":\n"
"    return ( 0.0";
    if( !(correctionFunctions.empty()) )
    {
      pythonCode << " + " << loopFactor  << " * ( ";
      for( std::vector< std::string >::iterator
           whichCorrection( correctionFunctions.begin() );
           correctionFunctions.end() > whichCorrection;
           ++whichCorrection )
        // each correction must be added:
      {
        if( whichCorrection > correctionFunctions.begin() )
        {
          pythonCode << " + ";
        }
        pythonCode << *whichCorrection << argumentsString;
      }
      pythonCode << " )";
    }
    pythonCode << " )\n"
"\n"
"def " << loopCorrectedPotentialFunctionName << argumentsString << ":\n"
"    return ( " << polynomialPartFunctionName << argumentsString << " + "
<< loopCorrectionFunctionName << argumentsString << " )";
  }

} /* namespace Vevacious */
