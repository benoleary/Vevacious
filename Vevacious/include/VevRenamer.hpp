/*
 * VevRenamer.hpp
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

#ifndef VEVRENAMER_HPP_
#define VEVRENAMER_HPP_

#include <string>
#include <map>
#include <vector>
#include <list>
#include <stdexcept>
#include "BOLlib/include/ArgumentParser.hpp"
#include "BOLlib/include/StringParser.hpp"

namespace Vevacious
{
  // this does the job of turning user-defined names for VEVs provided through
  // SARAH into single-character names as required by PyMinuit, & deals with
  // ordering issues.
  class VevRenamer
  {
  public:
    static std::string const emptyString;

    VevRenamer();
    ~VevRenamer();

    void setVevNames( std::map< std::string, std::string > const& inputVevs );
    void sortVevsToMatchSolutions(
                           std::vector< char > const& vevsOrderedBySolutions );
    std::string replaceUserVevNames( std::string const& inputString ) const;
    std::string const& getHumanReadableVevNameMap() const;
    char getInternalVevName( std::string const& userVevName ) const;
    std::string const& getInternalVevNamesAsQuotedCharList() const;
    std::string const& getInternalVevNamesAsUnquotedCharList() const;
    std::vector< char > const& getVevsOrderedBySolutions() const;
    std::pair< std::string, std::string >
    getUserVevNameWithValue( char const internalVevName ) const;
    std::string getUserVevName( char const internalVevName ) const
    { return getUserVevNameWithValue( internalVevName ).first; }
    std::string getInputVevValue( char const internalVevName ) const
    { return getUserVevNameWithValue( internalVevName ).second; }


  protected:
    static std::string const possibleInternalVevNames;

    static bool
    isOrderedLongToShort( std::pair< std::string, char > const& firstPair,
                          std::pair< std::string, char > const& secondPair );

    std::map< char, std::pair< std::string, std::string > >
    internalVevNamesToInputVevNamesAndValues;
    std::list< std::pair< std::string, char > >
    inputVevNamesToInternalVevNames;
    std::vector< char > internalVevNamesOrderedForSolutions;
    std::string humanReadableVevNameMap;
    std::string internalVevNamesAsQuotedCharList;
    std::string internalVevNamesAsUnquotedCharList;
    std::string currentVevName;
  };




  inline void VevRenamer::sortVevsToMatchSolutions(
                            std::vector< char > const& vevsOrderedBySolutions )
  {
    internalVevNamesOrderedForSolutions = vevsOrderedBySolutions;
  }

  inline std::string const& VevRenamer::getHumanReadableVevNameMap() const
  {
    return humanReadableVevNameMap;
  }

  inline char
  VevRenamer::getInternalVevName( std::string const& userVevName ) const
  {
    for( std::list< std::pair< std::string, char > >::const_iterator
         whichPair( inputVevNamesToInternalVevNames.begin() );
         inputVevNamesToInternalVevNames.end() != whichPair;
         ++whichPair )
    {
      if( 0 == userVevName.compare( whichPair->first ) )
      {
        return whichPair->second;
      }
    }
    return '?';
  }

  inline std::string const&
  VevRenamer::getInternalVevNamesAsQuotedCharList() const
  {
    return internalVevNamesAsQuotedCharList;
  }

  inline std::string const&
  VevRenamer::getInternalVevNamesAsUnquotedCharList() const
  {
    return internalVevNamesAsUnquotedCharList;
  }

  inline std::vector< char > const&
  VevRenamer::getVevsOrderedBySolutions() const
  {
    return internalVevNamesOrderedForSolutions;
  }

  inline bool VevRenamer::isOrderedLongToShort(
                               std::pair< std::string, char > const& firstPair,
                             std::pair< std::string, char > const& secondPair )
  {
    return ( firstPair.first.size() >= secondPair.first.size() );
  }

  inline std::pair< std::string, std::string >
  VevRenamer::getUserVevNameWithValue( char const internalVevName ) const
  {
    std::map< char, std::pair< std::string, std::string > >::const_iterator
    vevFinder( internalVevNamesToInputVevNamesAndValues.find(
                                                           internalVevName ) );
    if( internalVevNamesToInputVevNamesAndValues.end() == vevFinder )
    {
      std::string exceptionString( "No user-definition for \'" );
      exceptionString.append( 1,
                              internalVevName );
      exceptionString.append( "\'!" );
      throw std::out_of_range( exceptionString );
    }
    return vevFinder->second;
  }

} /* namespace Vevacious */
#endif /* VEVRENAMER_HPP_ */
