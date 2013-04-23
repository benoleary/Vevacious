/*
 * SarahInterpreter.hpp
 *
 *  Created on: Feb 7, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2013 Ben O'Leary
 *
 *      This file is part of Vevacious, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.Vevacious.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef SARAHINTERPRETER_HPP_
#define SARAHINTERPRETER_HPP_

#include <string>
#include <map>
#include "BOLlib/include/AsciiXmlParser.hpp"
#include "BOLlib/include/StringParser.hpp"
#include "BOLlib/include/VectorlikeArray.hpp"
#include "SarahSlhaConverter.hpp"
#include "VevRenamer.hpp"

namespace Vevacious
{
  // this class bundles together SarahSlhaConverter & VevRenamer for
  // convenience.
  class SarahInterpreter
  {
  public:
    static void removeNewlinesFrom( std::string& stringToEdit );

    SarahInterpreter();
    ~SarahInterpreter();

    void setSlhaFile( std::string const& slhaFilename );
    void setBlockPrefixes(
                   std::map< std::string, std::string > const& blockPrefixes );
    // this first sets the block prefixes to empty strings, *then* sets them to
    // the attributes given, so if "tree" is given, but "loop" is not,
    // treeBlockPrefix will be set to the given string, while loopBlockPrefix
    // will be set to an empty string.
    void setPositiveVevs( std::string const& positiveVevs );
    void setVevNames( std::map< std::string, std::string > const& inputVevs );
    void setVevScalingString( std::string const& vevScalingString );
    std::string const& operator()( std::string const& inputString,
                                   bool const useTreeRatherThanLoop = true,
                                   bool const writingTadpole = false );
    // this returns a reference to a string where all SLHA block elements in
    // inputString have been replaced by their values and all the user-defined
    // VEV names in inputString have been replaced by their internal
    // single-char names. if emptyIfHasSlhaZero is true and any of the SLHA
    // substitutions was for zero, a reference to an empty string is returned
    // instead. if emptyIfHasImaginaryUnit is true and any chars from
    // charsDenotingImaginaryUnit were in the string after accounting for SLHA
    // substitutions, a reference to an empty string is returned instead.
    // useTreeRatherThanLoop determines which prefix is prepended to SLHA block
    // names as a 1st replacement attempt, before the block name without a
    // prefix is used if the prepended block name is not found for the
    // appropriate index.
    void sortVevsToMatchSolutions(
                           std::vector< char > const& vevsOrderedBySolutions );
    double getSlhaScale() const;
    std::string const& getHumanReadableVevNameMap() const;
    std::string const& getInternalVevNamesAsQuotedCharList() const;
    std::string const& getInternalVevNamesAsUnquotedCharList() const;
    std::vector< char > const& getVevsOrderedBySolutions() const;
    std::string getUserVevName( char const internalVevName ) const;
    std::string getInputVevValue( char const internalVevName );
    std::string getInputVevPoint();
    std::string const& getPositiveInternalVevs() const;


  protected:
    static std::string const imaginaryUnitChars;

    SarahSlhaConverter* slhaConverter;
    std::string treeBlockPrefix;
    std::string loopBlockPrefix;
    VevRenamer vevRenamer;
    std::string vevScalingString;
    std::string argumentString;
    std::string returnString;
    size_t nextCharPosition;
    size_t lastCharPosition;
    std::string positiveInternalVevs;
  };





  inline void SarahInterpreter::removeNewlinesFrom( std::string& stringToEdit )
  {
    size_t newlinePosition( stringToEdit.find( '\n' ) );
    while( std::string::npos != newlinePosition )
    {
      stringToEdit.erase( newlinePosition,
                          1 );
      newlinePosition = stringToEdit.find( '\n',
                                           newlinePosition );
    }
  }

  inline void SarahInterpreter::setSlhaFile( std::string const& slhaFilename )
  {
    delete slhaConverter;
    slhaConverter = new SarahSlhaConverter( slhaFilename );
  }

  inline void SarahInterpreter::setBlockPrefixes(
                    std::map< std::string, std::string > const& blockPrefixes )
  // this first sets the block prefixes to empty strings, *then* sets them to
  // the attributes given, so if "tree" is given, but "loop" is not,
  // treeBlockPrefix will be set to the given string, while loopBlockPrefix
  // will be set to an empty string.
  {
    treeBlockPrefix.assign( "" );
    loopBlockPrefix.assign( "" );
    std::map< std::string, std::string >::const_iterator
    prefixFinder( blockPrefixes.find( "tree" ) );
    if( blockPrefixes.end() != prefixFinder )
    {
      treeBlockPrefix.assign( prefixFinder->second );
    }
    prefixFinder = blockPrefixes.find( "loop" );
    if( blockPrefixes.end() != prefixFinder )
    {
      loopBlockPrefix.assign( prefixFinder->second );
    }
  }


  inline void SarahInterpreter::setVevNames(
                        std::map< std::string, std::string > const& inputVevs )
  {
    vevRenamer.setVevNames( inputVevs );
  }

  inline void
  SarahInterpreter::setVevScalingString( std::string const& vevScalingString )
  {
    this->vevScalingString.assign( vevScalingString );
  }

  inline void SarahInterpreter::sortVevsToMatchSolutions(
                            std::vector< char > const& vevsOrderedBySolutions )
  {
    vevRenamer.sortVevsToMatchSolutions( vevsOrderedBySolutions );
  }

  inline double SarahInterpreter::getSlhaScale() const
  {
    return slhaConverter->getScale();
  }

  inline std::string const&
  SarahInterpreter::getHumanReadableVevNameMap() const
  {
    return vevRenamer.getHumanReadableVevNameMap();
  }

  inline std::string const&
  SarahInterpreter::getInternalVevNamesAsQuotedCharList() const
  {
    return vevRenamer.getInternalVevNamesAsQuotedCharList();
  }

  inline std::string const&
  SarahInterpreter::getInternalVevNamesAsUnquotedCharList() const
  {
    return vevRenamer.getInternalVevNamesAsUnquotedCharList();
  }

  inline std::vector< char > const&
  SarahInterpreter::getVevsOrderedBySolutions() const
  {
    return vevRenamer.getVevsOrderedBySolutions();
  }

  inline std::string
  SarahInterpreter::getUserVevName( char const internalVevName ) const
  {
    return vevRenamer.getUserVevName( internalVevName );
  }

  inline std::string
  SarahInterpreter::getInputVevValue( char const internalVevName )
  {
    return BOL::StringParser::trimFromFrontAndBack( (*slhaConverter)(
                              vevRenamer.getInputVevValue( internalVevName ) ),
                                                    "\"\'" );
  }

  inline std::string
  SarahInterpreter::getInputVevPoint()
  {
    std::string vevPoint( "" );
    for( std::vector< char >::const_iterator
         whichField( vevRenamer.getVevsOrderedBySolutions().begin() );
         vevRenamer.getVevsOrderedBySolutions().end() > whichField;
         ++whichField )
    {
      if( !(vevPoint.empty()) )
      {
        vevPoint.append( ", " );
      }
      vevPoint.append( BOL::StringParser::trimFromFrontAndBack(
                (*slhaConverter)( vevRenamer.getInputVevValue( *whichField ) ),
                                                                "\"\'" ) );
    }
    return vevPoint;
  }

  inline std::string const&
  SarahInterpreter::getPositiveInternalVevs() const
  {
    return positiveInternalVevs;
  }

} /* namespace Vevacious */
#endif /* SARAHINTERPRETER_HPP_ */
