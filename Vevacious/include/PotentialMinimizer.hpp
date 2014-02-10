/*
 * PotentialMinimizer.hpp
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

#ifndef POTENTIALMINIMIZER_HPP_
#define POTENTIALMINIMIZER_HPP_

#include <set>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include "BOLlib/include/ArgumentParser.hpp"
#include "BOLlib/include/AsciiXmlParser.hpp"
#include "BOLlib/include/StringParser.hpp"
#include "SLHA.hpp"
#include "SarahInterpreter.hpp"

namespace Vevacious
{
  // this class writes a Python program to use PyMinuit to minimize a
  // loop-corrected potential energy function from the form given by SARAH
  // combined with an SLHA file.
  class PotentialMinimizer
  {
  public:
    static std::string const namesOfVevs;
    static std::string const internalVevNamesToUserNames;
    static std::string const vevsTakenPositive;
    static std::string const energyScale;
    static std::string const energyScaleFourth;
    static std::string const userVevsAsMathematica;
    static std::string const userVevsAsXml;
    static std::string const vevDictionaryToArray;
    static std::string const vevOrigin;
    static std::string const inputVevsPoint;
    static std::string const functionFromDictionary;
    static std::string const functionFromArray;
    static std::string const loopCorrectedPotential;
    static std::string const treeLevelPotential;

    PotentialMinimizer();
    ~PotentialMinimizer();

    void setPolynomialPartOfPotential(
                                std::string const& polynomialPartOfPotential );
    void addMassSquaredMatrix(
                      std::map< std::string, std::string > const& attributeMap,
                               std::string const& massSquaredMatrix );
    std::string
    prepareParameterDependentPython( SarahInterpreter& sarahInterpreter );


  protected:
    typedef std::pair< std::string, std::string > StringPair;
    typedef std::pair< StringPair, StringPair > StringQuadruple;

    static std::string const overallFactorAttributeName;
    static std::string const subtractionConstantAttributeName;
    static std::string const spinTypeAttributeName;
    static std::string const defaultSubtractionConstantString;

    static bool isOrderedLongToShort(
                        std::pair< std::string, std::string > const& firstPair,
                     std::pair< std::string, std::string > const& secondPair );

    BOL::AsciiXmlParser sarahParser;
    std::map< std::string, std::string >::const_iterator attributeFinder;
    std::stringstream pythonCode;
    std::string polynomialPartOfPotential;
    std::vector< StringQuadruple > factoredMassSquaredMatrices;
    std::vector< std::string > massesSquaredFunctions;
    std::stringstream stringParser;
    std::string argumentsWithoutBrackets;
    std::string argumentsWithBrackets;

    void preparePolynomialPartOfPotential( SarahInterpreter& sarahInterpreter );
    void prepareLoopCorrections( SarahInterpreter& sarahInterpreter );
    void prepareMinimizations( SarahInterpreter& sarahInterpreter );
  };





  inline void PotentialMinimizer::setPolynomialPartOfPotential(
                                 std::string const& polynomialPartOfPotential )
  {
    this->polynomialPartOfPotential.assign( polynomialPartOfPotential );
    SarahInterpreter::removeNewlinesFrom( this->polynomialPartOfPotential );
  }

  inline std::string
  PotentialMinimizer::prepareParameterDependentPython(
                                           SarahInterpreter& sarahInterpreter )
  {
    // the function names should explain the process well enough...
    argumentsWithoutBrackets.assign(
                    sarahInterpreter.getInternalVevNamesAsUnquotedCharList() );
    argumentsWithBrackets.assign( "( " + argumentsWithoutBrackets + " )" );
    preparePolynomialPartOfPotential( sarahInterpreter );
    prepareLoopCorrections( sarahInterpreter );
    return pythonCode.str();
  }

  inline void PotentialMinimizer::preparePolynomialPartOfPotential(
                                           SarahInterpreter& sarahInterpreter )
  {
    if( polynomialPartOfPotential.empty() )
    {
      throw std::runtime_error( "no polynomial part of potential" );
    }
    // 2 temporary strings to avoid problems with function evaluation timings.
    std::string
    specificTreeLevelPotential( sarahInterpreter( polynomialPartOfPotential,
                                                  true ) );
    std::string
    specificPolynomialPart( sarahInterpreter( polynomialPartOfPotential,
                                              false ) );
    pythonCode
    << "def TreeLevelPotential( " << argumentsWithoutBrackets
    << ", temperatureValue = 0.0 ):\n"
    << "    return ( " << specificTreeLevelPotential << " )\n"
    << "\n"
    << "def PolynomimalPartOfPotential" << argumentsWithBrackets << ":\n"
    << "    return ( " << specificPolynomialPart << " )\n"
    << "\n";
  }

} /* namespace Vevacious */
#endif /* POTENTIALMINIMIZER_HPP_ */
