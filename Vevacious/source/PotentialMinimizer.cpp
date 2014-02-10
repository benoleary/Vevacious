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
  std::string const PotentialMinimizer::overallFactorAttributeName( "factor" );
  std::string const
  PotentialMinimizer::subtractionConstantAttributeName( "constant" );
  std::string const
  PotentialMinimizer::spinTypeAttributeName( "spin" );
  std::string const
  PotentialMinimizer::defaultSubtractionConstantString( "1.5" );


  PotentialMinimizer::PotentialMinimizer() :
    sarahParser( false ),
    attributeFinder(),
    pythonCode(),
    polynomialPartOfPotential( "" ),
    factoredMassSquaredMatrices(),
    massesSquaredFunctions(),
    stringParser(),
    argumentsWithoutBrackets( "" ),
    argumentsWithBrackets( "()" )
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
                      std::map< std::string, std::string > const& attributeMap,
                                         std::string const& massSquaredMatrix )
  {
    attributeFinder
    = attributeMap.find( subtractionConstantAttributeName );
    std::string subtractionConstant( defaultSubtractionConstantString );
    if( attributeMap.end() != attributeFinder )
    {
      subtractionConstant.assign( BOL::StringParser::trimFromFrontAndBack(
                                                       attributeFinder->second,
                                                                    "\"\'" ) );
    }
    attributeFinder = attributeMap.find( spinTypeAttributeName );
    std::string spinType( "" );
    if( attributeMap.end() != attributeFinder )
    {
      spinType.assign( BOL::StringParser::trimFromFrontAndBack(
                                                       attributeFinder->second,
                                                                    "\"\'" ) );
    }
    attributeFinder = attributeMap.find( overallFactorAttributeName );
    if( attributeMap.end() != attributeFinder )
    {
      factoredMassSquaredMatrices.push_back( StringQuadruple() );
      factoredMassSquaredMatrices.back().first.first.assign(
              BOL::StringParser::trimFromFrontAndBack( attributeFinder->second,
                                                       "\"\'" ) );
      factoredMassSquaredMatrices.back().first.second.assign(
                                                         subtractionConstant );
      factoredMassSquaredMatrices.back().second.first.assign( spinType );
      factoredMassSquaredMatrices.back().second.second.assign(
                                                           massSquaredMatrix );
      SarahInterpreter::removeNewlinesFrom(
                            factoredMassSquaredMatrices.back().second.second );
    }
  }

  void PotentialMinimizer::prepareLoopCorrections(
                                           SarahInterpreter& sarahInterpreter )
  {
    massesSquaredFunctions.clear();
    std::string massSquaredMatrix( "" );
    std::string overallFactor( "" );
    std::string subtractionConstant( "" );
    std::string massSquaredFunctionName( "" );
    BOL::VectorlikeArray< std::string > massTerms;
    for( std::vector< StringQuadruple >::iterator
         whichQuadruple( factoredMassSquaredMatrices.begin() );
         factoredMassSquaredMatrices.end() > whichQuadruple;
         ++whichQuadruple )
    {
      massTerms.clearEntries();
      BOL::StringParser::parseByChar( sarahInterpreter(
                                               whichQuadruple->second.second ),
                                      massTerms,
                                      ',' );
      stringParser.clear();
      stringParser.str( "" );
      stringParser
      << "MassesSquared" << ( massesSquaredFunctions.size() + 1 );
      massSquaredFunctionName.assign( stringParser.str() );
      pythonCode <<
"# The masses-squared from this function are in GeV^2:\n"
"def " << massSquaredFunctionName << argumentsWithBrackets << ":\n"
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
"    massSquaredMatrix = massSquaredElementArray.reshape( math.sqrt(\n"
"                                       massSquaredElementArray.size ), -1 )\n"
"    return numpy.linalg.eigvalsh( massSquaredMatrix )\n"
"\n";
      stringParser
      << ", " << whichQuadruple->first.first
      << ", " << whichQuadruple->first.second
      << ", " << whichQuadruple->second.first << " ]";
      massesSquaredFunctions.push_back( "[ " + stringParser.str() );
      // all the corrections have to be noted for later.
    }
    pythonCode <<
"\n"
"\n"
"MassesSquareds = [ \n";
    for( std::vector< std::string >::iterator
         whichMassesSquared( massesSquaredFunctions.begin() );
         massesSquaredFunctions.end() > whichMassesSquared;
         ++whichMassesSquared )
    {
      if( whichMassesSquared > massesSquaredFunctions.begin() )
      {
        pythonCode << ",\n                   ";
      }
      pythonCode << *whichMassesSquared;
    }
    pythonCode << " ]\n"
"\n"
"def LoopCorrectedPotential( " << argumentsWithoutBrackets
<<                                              ", temperatureValue = 0.0 ):\n"
"    loopCorrections = 0.0\n"
"    for MassesSquared in MassesSquareds:\n"
"        loopCorrections += MassSquaredCorrections( MassesSquared[ 0 ]"
<<                                               argumentsWithBrackets << ",\n"
"                                                   MassesSquared[ 1 ],\n"
"                                                   MassesSquared[ 2 ] )\n"
"    return ( PolynomimalPartOfPotential" << argumentsWithBrackets << "\n"
"             + ( loopFactor * loopCorrections ) )\n"
"\n"
"\n"
"def LoopAndThermalCorrectedPotential( " << argumentsWithoutBrackets
<<                                                    ", temperatureValue ):\n"
"    loopCorrections = 0.0\n"
"    thermalCorrections = 0.0\n"
"    temperatureSquared = ( temperatureValue * temperatureValue )\n"
"    for MassesSquared in MassesSquareds:\n"
"        massesSquaredValues = MassesSquared[ 0 ]"
<<                                                argumentsWithBrackets << "\n"
"        loopCorrections += MassSquaredCorrections( massesSquaredValues,\n"
"                                                   MassesSquared[ 1 ],\n"
"                                                   MassesSquared[ 2 ] )\n"
"        if ( temperatureValue > 0.0 ):\n"
"            adjustedOverallFactor = MassesSquared[ 1 ]\n"
"            if ( MassesSquared[ 3 ] is \"vector\" ):\n"
"                adjustedOverallFactor = ( ( 2.0 * adjustedOverallFactor )\n"
"                                          / 3.0 )\n"
"            thermalCorrections += ThermalCorrections( massesSquaredValues,\n"
"                                                     adjustedOverallFactor,\n"
"                                                      temperatureSquared )\n"
"    return ( PolynomimalPartOfPotential" << argumentsWithBrackets << "\n"
"             + ( loopFactor * loopCorrections )\n"
"          + ( thermalFactor * temperatureValue**4 * thermalCorrections ) )\n"
"\n\n"
"# PotentialScaler class:\n"
"\n"
"class PotentialScaler:\n"
"    \"\"\"\n"
"    This class exists solely to provide a function that is a scaled version\n"
"    of the function given to the constructor for PyMinuit, and a set of\n"
"    functions using arrays that the CosmoTransitions objects want to use.\n"
"    \"\"\"\n"
"\n"
"    def __init__( self, FunctionToScale, temperatureValue = 0.0 ):\n"
"        self.FunctionToScale = FunctionToScale\n"
"        self.functionAtOrigin = FunctionToScale( ";
    std::vector< char > const&
    vevsOrderedBySolution( sarahInterpreter.getVevsOrderedBySolutions() );
    for( unsigned int whichVev( 0 );
        vevsOrderedBySolution.size() > whichVev;
        ++whichVev )
    {
      if( 0 < whichVev )
      {
        pythonCode << ", ";
      }
      pythonCode << "0.0";
    }
    pythonCode << " )\n"
"        self.temperatureValue = temperatureValue\n"
"        \n"
"    def ScaledFunctionFromScaledArguments( self, " << argumentsWithoutBrackets
<<                                                                      " ):\n"
"        functionValue = ( self.FunctionToScale( ";
    for( std::vector< char >::const_iterator
         whichName( vevsOrderedBySolution.begin() );
         vevsOrderedBySolution.end() > whichName;
         ++whichName )
    {
      pythonCode << "( " << *whichName << " * energyScale )\n"
      << "                                                ";
    }
    pythonCode << "self.temperatureValue )\n"
"                          - self.functionAtOrigin )\n"
"        return ( inverseScaleFourthed\n"
"                 * functionValue )\n"
"\n"
"    def PotentialFromArray( self, pointAsArray ):\n"
"        return FunctionFromArray( self.FunctionToScale,\n"
"                                  pointAsArray,\n"
"                                  self.temperatureValue )\n"
"\n"
"    def PotentialFromMatrix( self, arrayOfArrays ):\n"
"        if ( ( numberOfFields, ) == arrayOfArrays.shape ):\n"
"            return self.PotentialFromArray( arrayOfArrays )\n"
"        elif ( ( len( arrayOfArrays ), numberOfFields )\n"
"               == arrayOfArrays.shape ):\n"
"            returnArray = numpy.zeros( len( arrayOfArrays ) )\n"
"            for whichIndex in range( len( arrayOfArrays ) ):\n"
"                returnArray[ whichIndex ] = self.PotentialFromArray(\n"
"                                              arrayOfArrays[ whichIndex ] )\n"
"            return returnArray\n"
"        else:\n"
"            return None\n"
"\n"
"    def GradientFromArray( self, pointAsArray ):\n"
"        potentialAtPoint = self.PotentialFromArray( pointAsArray )\n"
"        gradientArray = numpy.zeros( len( pointAsArray ) )\n"
"        for whichField in range( len( pointAsArray ) ):\n"
"            displacedPoint = pointAsArray.copy()\n"
"            displacedPoint[ whichField ] += numericalStepSize\n"
"            gradientArray[ whichField ] = ( ( self.PotentialFromArray(\n"
"                                                           displacedPoint )\n"
"                                              - potentialAtPoint )\n"
"                                            / numericalStepSize )\n"
"        return gradientArray\n"
"\n"
"    def GradientFromMatrix( self, arrayOfArrays ):\n"
"        if ( ( numberOfFields, ) == arrayOfArrays.shape ):\n"
"            return self.GradientFromArray( arrayOfArrays )\n"
"        elif ( ( len( arrayOfArrays ), numberOfFields )\n"
"               == arrayOfArrays.shape ):\n"
"            returnMatrix = arrayOfArrays.copy()\n"
"            for whichIndex in range( len( arrayOfArrays ) ):\n"
"                returnMatrix[ whichIndex ] = self.GradientFromArray(\n"
"                                              arrayOfArrays[ whichIndex ] )\n"
"            return returnMatrix\n"
"        else:\n"
"            return None\n"
"\n"
"# End of PotentialScaler class\n"
"\n";
  }

} /* namespace Vevacious */
