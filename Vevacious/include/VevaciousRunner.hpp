/*
 * VevaciousRunner.hpp
 *
 *  Created on: Mar 5, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2013 Ben O'Leary
 *
 *      This file is part of Vevacious, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.Vevacious.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef VEVACIOUSRUNNER_HPP_
#define VEVACIOUSRUNNER_HPP_

#include <fstream>
#include <stdexcept>
#include <sstream>
#include "BOLlib/include/ArgumentParser.hpp"
#include "BOLlib/include/AsciiXmlParser.hpp"
#include "BOLlib/include/UsefulStuff.hpp"
#include "BOLlib/include/WaitingOnSubprocessExecutor.hpp"
#include "SarahInterpreter.hpp"
#include "TadpoleSolver.hpp"
#include "PotentialMinimizer.hpp"

namespace Vevacious
{
  // this class wraps all the Vevacious functionality together.
  class VevaciousRunner
  {
  public:
    static std::string const vevaciousVersion;
    static std::string const vevaciousDocumentation;

    VevaciousRunner( std::string const& modelFilename,
                     std::string const& hom4ps2Directory,
                     std::string const& homotopyType,
                     std::string const& pythonMinuitWrapper );
    VevaciousRunner( BOL::ArgumentParser& argumentParser );
    ~VevaciousRunner();

    void setResultsFilename( std::string const& resultsFilename )
    { this->resultsFilename.assign( resultsFilename ); }
    void setImaginaryPartTolerance( double imaginaryTolerance )
    { this->imaginaryTolerance = imaginaryTolerance; }
    void setImaginaryPartTolerance( std::string const& imaginaryTolerance )
    { this->imaginaryTolerance
      = BOL::StringParser::stringToDouble( imaginaryTolerance ); }
    void
    setMinuitNudgesOffSaddlePoints( std::vector< double > saddleNudgeList )
    { this->saddleNudgeList = saddleNudgeList; }
    void setMinuitNudgesOffSaddlePoints(
                               std::string const& nudgesAsCommaSeparatedList );
    void setMaximumMinuitNudgesOffSaddlePoints( int maximumSaddleSplits )
    { saddleNudgeList.resize( maximumSaddleSplits,
                              saddleNudgeList.back() ); }
    void setMaximumMinuitNudgesOffSaddlePoints(
                                       std::string const& maximumSaddleSplits )
    { setMaximumMinuitNudgesOffSaddlePoints(
                     BOL::StringParser::stringToInt( maximumSaddleSplits ) ); }
    void setMinuitRollingTolerance( double rollingTolerance )
    { this->rollingTolerance = rollingTolerance; }
    void setMinuitRollingTolerance( std::string const& rollingTolerance )
    { this->rollingTolerance
      = BOL::StringParser::stringToDouble( rollingTolerance ); }
    void setPathToCosmotransitions( std::string const& pathToCosmotransitions )
    { this->pathToCosmotransitions.assign( pathToCosmotransitions ); }
    void setLifetimeThreshold( double lifetimeThreshold )
    { this->lifetimeThreshold = lifetimeThreshold; }
    void setLifetimeThreshold( std::string const& lifetimeThreshold )
    { setLifetimeThreshold(
                    BOL::StringParser::stringToDouble( lifetimeThreshold ) ); }
    void setThermalSurvivalThreshold( double survivalThreshold )
    { this->survivalThreshold = survivalThreshold; }
    void setThermalSurvivalThreshold( std::string const& survivalThreshold )
    { setThermalSurvivalThreshold(
                    BOL::StringParser::stringToDouble( survivalThreshold ) ); }
    void setShouldTunnel( bool const shouldTunnel )
    { if( shouldTunnel ){ this->shouldTunnel = "True"; }
      else{ this->shouldTunnel = "False"; } }
    void setShouldTunnelThermally( bool const tunnelThermally )
    { if( tunnelThermally ){ this->tunnelThermally = "True"; }
      else{ this->tunnelThermally = "False"; } }
    void setShouldDeformTunnelPaths( bool const deformTunnelPaths )
    { if( deformTunnelPaths ){ this->deformTunnelPaths = "True"; }
      else{ this->deformTunnelPaths = "False"; } }
    void appendTreeLevelExtrema( std::string const& slhaFilename,
                                 std::string const solutionsFilename
                                 = "VevaciousTreeLevelExtrema.py" );
    void findTreeLevelExtrema( std::string const& slhaFilename,
                               std::string const solutionsFilename
                               = "VevaciousTreeLevelExtrema.py" )
    { overwriteTreeLevelExtrema( slhaFilename,
                                 solutionsFilename ); }
    void overwriteTreeLevelExtrema( std::string const& slhaFilename,
                                    std::string const solutionsFilename
                                    = "VevaciousTreeLevelExtrema.py" );
    void writeCalculatedTreeLevelExtrema( std::string const& outputFilename );
    void prepareParameterDependentPython( std::string const& slhaFilename,
                                          std::string const pythonFilename
                                          = "VevaciousParameterDependent.py" );
    void runPython( std::string const pythonFilename = "Vevacious.py" );
    // this checks to see if there is a file called pythonFilename, then, if it
    // exists, runs it with the terminal's "python". if it doesn't exist, then
    // the default file with name given by defaultPythonFilename is written &
    // run.
    bool runPython( std::string const pythonFilename,
                    int const deformationPatience );
    // this does the same as runPython( std::string const pythonFilename ), but
    // killing the Python subprocess after deformationPatience seconds if it
    // takes that long. it returns true if the subprocess finished without
    // having to be killed, false if the subprocess took too long.
    void appendResultsToSlha( std::string const& slhaFilename,
                   bool const properFormatRatherThanSspReadable = true ) const;
    // this takes the XML-format results from the file called resultsFilename,
    // & appends the results & warnings in custom SLHA blocks to the end of the
    // file name slhaFilename.


  protected:
    static std::string const defaultPythonFilename;
    static std::string const treeLevelExtremaHeader;
    static BOL::StringParser const slhaIndexMaker;
    static BOL::StringParser const slhaDoubleMaker;

    static double doubleFromQuotedString( std::string const& quotedString );
    static std::string writeVevaciousResultsBlockLine( int const firstIndex,
                                                       int const secondIndex,
                                                       double floatValue,
                                                std::string const& givenString,
                                bool const properFormatRatherThanSspReadable );

    SarahInterpreter sarahInterpreter;
    TadpoleSolver tadpoleSolver;
    PotentialMinimizer potentialMinimizer;
    std::string pythonMinuitWrapper;
    std::stringstream treeLevelExtrema;
    bool firstWriteOfExtrema;
    std::string resultsFilename;
    double imaginaryTolerance;
    std::vector< double > saddleNudgeList;
    double rollingTolerance;
    std::string pathToCosmotransitions;
    double lifetimeThreshold;
    double survivalThreshold;
    std::string shouldTunnel;
    std::string tunnelThermally;
    std::string deformTunnelPaths;

    void calculateTreeLevelExtrema( std::string const& slhaFilename );
    void writeDefaultPythonProgram() const;
    std::string
    prepareToRunPythonCommand( std::string const& pythonFilename ) const;
  };





  inline void VevaciousRunner::setMinuitNudgesOffSaddlePoints(
                                std::string const& nudgesAsCommaSeparatedList )
  {
    saddleNudgeList.clear();
    std::string spaceSeparatedString( BOL::StringParser::trimFromFrontAndBack(
                                                    nudgesAsCommaSeparatedList,
                    ( BOL::StringParser::whitespaceAndNewlineChars + "," ) ) );
    BOL::StringParser::substituteCharacterWith( spaceSeparatedString,
                                                ',',
                                                ' ' );
    if( !(spaceSeparatedString.empty()) )
    {
      std::stringstream streamToParse( spaceSeparatedString );
      double parsedIntAsDouble;
      while( streamToParse.good() )
      {
        streamToParse >> parsedIntAsDouble;
        saddleNudgeList.push_back( parsedIntAsDouble );
      }
    }
  }

  inline void
  VevaciousRunner::appendTreeLevelExtrema( std::string const& slhaFilename,
                                          std::string const solutionsFilename )
  {
    BOL::UsefulStuff::runSystemCommand( "rm " + solutionsFilename + "c" );
    calculateTreeLevelExtrema( slhaFilename );
    std::ofstream outputFile( solutionsFilename.c_str(),
                              std::ios::app );
    if( !(outputFile.good()) )
    {
      throw std::runtime_error( "could not open \"" + solutionsFilename
                                + "\" to write tree-level extrema." );
    }
    if( firstWriteOfExtrema )
    {
      outputFile << treeLevelExtremaHeader;
      firstWriteOfExtrema = false;
    }
    outputFile << treeLevelExtrema.str();
    outputFile.close();
  }

  inline void
  VevaciousRunner::overwriteTreeLevelExtrema( std::string const& slhaFilename,
                                          std::string const solutionsFilename )
  {
    BOL::UsefulStuff::runSystemCommand( "rm " + solutionsFilename + "c" );
    calculateTreeLevelExtrema( slhaFilename );
    std::ofstream outputFile( solutionsFilename.c_str() );
    if( !(outputFile.good()) )
    {
      throw std::runtime_error( "could not open \"" + solutionsFilename
                                + "\" to write tree-level extrema." );
    }
    outputFile << treeLevelExtremaHeader << treeLevelExtrema.str();
    outputFile.close();
  }

  inline void VevaciousRunner::runPython( std::string const pythonFilename )
  // this checks to see if there is a file called pythonFilename, then, if it
  // exists, runs it with the terminal's "python". if it doesn't exist, then
  // the default file with name given by defaultPythonFilename is written &
  // run.
  {
    BOL::UsefulStuff::runSystemCommand(
                                 prepareToRunPythonCommand( pythonFilename ) );
  }

  inline bool VevaciousRunner::runPython( std::string const pythonFilename,
                                          int const deformationPatience )
   // this does the same as runPython( std::string const pythonFilename ), but
   // killing the Python subprocess after deformationPatience seconds if it
   // takes that long. it returns true if the subprocess finished without
   // having to be killed, false if the subprocess took too long.
   {
     std::string systemCommand( prepareToRunPythonCommand( pythonFilename ) );
     if( ( 0 < deformationPatience )
         &&
         ( 0.0 < lifetimeThreshold ) )
     {
       std::string forkedCommand( "bash " + systemCommand );
       BOL::WaitingOnSubprocessExecutor subprocessExecutor( forkedCommand,
                                                            true,
                                               ( 1000 * deformationPatience ) );
       // BOL::WaitingOnSubprocessExecutor takes a time in milliseconds, while
       // deformedPathPatience is in seconds.
       return subprocessExecutor.forkAndExecvAndWait();
     }
     BOL::UsefulStuff::runSystemCommand( systemCommand );
     return true;
   }

  inline double VevaciousRunner::doubleFromQuotedString(
                                              std::string const& quotedString )
  {
    return BOL::StringParser::stringToDouble(
                         BOL::StringParser::trimFromFrontAndBack( quotedString,
                                                                  "\'\"" ) );
  }

  inline std::string
  VevaciousRunner::writeVevaciousResultsBlockLine( int const firstIndex,
                                                   int const secondIndex,
                                                   double const floatValue,
                                         std::string const& givenString,
                                 bool const properFormatRatherThanSspReadable )
  {
    std::string returnString( "  " );
    returnString.append( slhaIndexMaker.intToString( firstIndex ) );
    returnString.append( "  " );
    returnString.append( slhaIndexMaker.intToString( secondIndex ) );
    returnString.append( "    " );
    returnString.append( slhaIndexMaker.doubleToString( floatValue ) );
    if( !properFormatRatherThanSspReadable )
    {
      returnString.append( " # " );
    }
    returnString.append( "    " );
    returnString.append( givenString );
    return returnString;
  }

  inline std::string VevaciousRunner::prepareToRunPythonCommand(
                                      std::string const& pythonFilename ) const
  {
    std::string systemCommand( "rm " + pythonFilename + "c" );
    BOL::UsefulStuff::runSystemCommand( systemCommand );
    systemCommand.assign( "python " );
    if( !(BOL::UsefulStuff::fileExists( pythonFilename )) )
    {
      writeDefaultPythonProgram();
      systemCommand.append( defaultPythonFilename );
    }
    else
    {
      systemCommand.append( pythonFilename );
    }
    return systemCommand;
  }

} /* namespace Vevacious */
#endif /* VEVACIOUSRUNNER_HPP_ */
