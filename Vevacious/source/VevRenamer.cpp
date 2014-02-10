/*
 * VevRenamer.cpp
 *
 *  Created on: Feb 7, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of Vevacious, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.Vevacious.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "../include/VevRenamer.hpp"

namespace Vevacious
{
  std::string const VevRenamer::emptyString( "" );
  std::string const VevRenamer::possibleInternalVevNames(
                              "abcdfghklmnpqrstuvwxyzABCDFGHKLMNPQRSTUVWXYZ" );

  VevRenamer::VevRenamer() :
    internalVevNamesToInputVevNamesAndValues(),
    inputVevNamesToInternalVevNames(),
    internalVevNamesOrderedForSolutions(),
    humanReadableVevNameMap( "" ),
    internalVevNamesAsQuotedCharList( "" ),
    internalVevNamesAsUnquotedCharList( "" ),
    currentVevName( "" )
  {
    // just an initialization list.
  }

  VevRenamer::~VevRenamer()
  {
    // does nothing.
  }


  void VevRenamer::setVevNames(
                        std::map< std::string, std::string > const& inputVevs )
  {
    if( possibleInternalVevNames.size() < inputVevs.size() )
    {
      std::cout
      << std::endl
      << "Error! Due to some unfortunate function variable name restrictions,"
      << " Vevacious cannot deal with more than "
      << possibleInternalVevNames.size() << " possible non-zero VEVs. ("
      << inputVevs.size() << " VEVs were given.)";
      std::cout << std::endl;

      throw std::overflow_error( "too many VEVs" );
    }
    internalVevNamesOrderedForSolutions.clear();
    char const* internalVevName( &(possibleInternalVevNames[ 0 ]) );
    humanReadableVevNameMap.assign( "" );
    internalVevNamesAsQuotedCharList.assign( "" );
    internalVevNamesAsUnquotedCharList.assign( "" );
    for( std::map< std::string, std::string >::const_iterator
         whichVev( inputVevs.begin() );
         inputVevs.end() != whichVev;
         ++whichVev )
    {
      currentVevName.assign( BOL::StringParser::trimFromFrontAndBack(
                                                               whichVev->first,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
      if( 2 > currentVevName.size() )
      {
        std::cout
        << std::endl
        << "ERROR! Unfortunately Vevacious needs all VEVs to have names that"
        << " are at least 2 characters in length each.";
        std::cout << std::endl;

        throw std::underflow_error(
                          "VEV name too short (\"" + currentVevName + "\")" );
      }
      internalVevNamesToInputVevNamesAndValues.insert( std::make_pair(
                                                              *internalVevName,
                                                                std::make_pair(
                                                                currentVevName,
                                                        whichVev->second ) ) );
      inputVevNamesToInternalVevNames.push_back( std::make_pair(
                                                                currentVevName,
                                                          *internalVevName ) );
      if( inputVevs.begin() != whichVev)
      {
        humanReadableVevNameMap.append( ", " );
        internalVevNamesAsQuotedCharList.append( ", " );
        internalVevNamesAsUnquotedCharList.append( ", " );
      }
      humanReadableVevNameMap.append( whichVev->first + " -> " );
      humanReadableVevNameMap.append( 1,
                                      *internalVevName );
      internalVevNamesAsQuotedCharList.append( "\'" );
      internalVevNamesAsQuotedCharList.append( 1,
                                               *internalVevName );
      internalVevNamesAsQuotedCharList.append( "\'" );
      internalVevNamesAsUnquotedCharList.append( 1,
                                                 *internalVevName );
      internalVevNamesOrderedForSolutions.push_back( *internalVevName );
      internalVevName += 1;
    }
    inputVevNamesToInternalVevNames.sort( &isOrderedLongToShort );
  }

  std::string
  VevRenamer::replaceUserVevNames( std::string const& inputString ) const
  {
    std::string returnString( inputString );
    size_t vevNamePosition( 0 );
    for( std::list< std::pair< std::string, char > >::const_iterator
         whichVev( inputVevNamesToInternalVevNames.begin() );
         inputVevNamesToInternalVevNames.end() != whichVev;
         ++whichVev )
    {
      vevNamePosition = returnString.find( whichVev->first );
      while( std::string::npos != vevNamePosition )
      {
        returnString.replace( vevNamePosition,
                              whichVev->first.size(),
                              1,
                              whichVev->second );
        vevNamePosition = returnString.find( whichVev->first,
                                             vevNamePosition );
      }
    }
    return returnString;
  }

} /* namespace Vevacious */
