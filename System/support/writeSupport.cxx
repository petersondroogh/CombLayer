/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   support/writeSupport.cxx
 *
 * Copyright (c) 2004-2018 by Stuart Ansell
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 *
 ****************************************************************************/
#include <iostream>
#include <iterator>
#include <iomanip>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cctype>
#include <complex>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <functional>
#include <boost/format.hpp>

#include "Exception.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Quaternion.h"

#include "support.h"
#include "writeSupport.h"

/*! 
  \file writeSupport.cxx 
*/

namespace  StrFunc
{

std::string
flukaNumber(const double D)
  /*!
    Process a number into a fluka style string
    \param D :: Number to use
  */
{
  boost::format FMTnum("%1$10.5g");
  return (FMTnum % D).str();
}

void
writeFLUKA(const std::string& Line,std::ostream& OX)
  /*!
    Write out the line in the fixed FLUKA format WHAT(1-6).
    Replace " - " by space to write empty WHAT cards.
    \param Line :: full FLUKA line
    \param OX :: ostream to write to
  */
{
  // this is expensive
  std::istringstream iss(Line);
  std::vector<std::string> whats(std::istream_iterator<std::string>{iss},
				 std::istream_iterator<std::string>());

  size_t i(1);
  for (std::string& w : whats)
    {
      if (w=="-") w=" ";

      if (w.size()>10)
	{
	  double D;
	  if (StrFunc::convert(w,D))
	    w=flukaNumber(D);
	  else
	    throw ColErr::InvalidLine(w,"String to long for FLUKA");
	}
      if (i % 8==0)
	OX<<std::setw(10)<<std::left<<w<<std::endl;
      else if (i % 8==1)
	OX<<std::setw(10)<<std::left<<w;	
      else 
	OX<<std::setw(10)<<std::right<<w;	
      i++;
    }
  
  if (i % 8 != 1) OX<<std::endl;
  return;
}

void
writeFLUKAhead(const std::string& HeadUnit,
	       const std::string& EndUnit,
	       const std::string& Line,
	       std::ostream& OX)
  /*!
    Write out the line in the fixed FLUKA format WHAT(1-6).
    Replace " - " by space to write empty WHAT cards.
    \param HeadUnit :: head name
    \param EndUnit :: end name
    \param Line :: full FLUKA line
    \param OX :: ostream to write to
  */
{
  const size_t NItem(6);
  // this is expensive
  std::istringstream iss(Line);
  std::vector<std::string> whats(std::istream_iterator<std::string>{iss},
				 std::istream_iterator<std::string>());

  size_t i(0);
  for (const std::string& w : whats)
    {
      if (!i) OX<<std::setw(10)<<std::left<<HeadUnit;
      OX<<std::setw(10)<<std::right<<w;
      i++;
      if (i==NItem)
	{
	  OX<<EndUnit<<std::endl;
	  i=0;
	}
    }
  
  if (i)
    {
      for(;i<NItem;i++)
	OX<<std::string(10,' ');
      OX<<EndUnit<<std::endl;
    }
  
  return;
}


void
writeMCNPX(const std::string& Line,std::ostream& OX)
/*!
  Write out the line in the limited form for MCNPX
  ie initial line from 0::72 after that 8 to 72
  (split on a space or comma)
  \param Line :: full MCNPX line
  \param OX :: ostream to write to
*/
{
  writeControl(Line,OX,72,8);
  return;
}

void
writeMCNPXcont(const std::string& Line,std::ostream& OX)
/*!
  Write out the line in the limited form for MCNPX
  ie initial line from 0::72 after that 8 to 72
  as if it is a split line
  \param Line :: full MCNPX line
  \param OX :: ostream to write to
*/
{
  writeControl(Line,OX,72,-8);
  return;
}

void
writeControl(const std::string& Line,std::ostream& OX,
	     const size_t LNmax,int insertDepth)
/*!
  Write out the line in the limited form for MCNPX
  ie initial line from 0::72 after that 8 to 72
  (split on a space or comma)
  \param Line :: full MCNPX line
  \param OX :: ostream to write to
  \param LNmax :: Maximium char count in a line
  \param insertDepth :: second line insert depth [-ve to include first line]
*/
{
  // Line
  size_t spcLen(0);
  if (insertDepth<0)
    {
      insertDepth *= -1;
      spcLen=static_cast<size_t>(insertDepth);
    }

  std::string::size_type pos(0);
  std::string X=Line.substr(0,LNmax-spcLen);    
  std::string::size_type posB=X.find_last_of(" ,");
  while(X.length() == LNmax-spcLen &&
	posB!=std::string::npos)
    {
      pos+=posB+1;
      if (!isspace(X[posB])) posB++;  // skip pass comma 
      X=fullBlock(X.substr(0,posB));
      if (!isEmpty(X))
	OX<<std::string(spcLen,' ')<<X<<std::endl;

      spcLen=static_cast<size_t>(insertDepth);
      X=Line.substr(pos,LNmax-spcLen);
      posB=X.find_last_of(" ,");
    }
    
  X=fullBlock(X);
  if (!isEmpty(X))
    OX<<std::string(spcLen,' ')<<X<<std::endl;
  return;
}

void
writeMCNPXcomment(const std::string& Line,std::ostream& OX)
/*!
  Write out the line in the limited form for MCNPX
  ie initial line from 0->72 after that 8 to 72
  (split on a space or comma). It puts int all in a comment
  so lines are cut by 2
  \param Line :: full MCNPX line
  \param OX :: ostream to write to
*/
{
  const size_t MaxLine(72);

  std::string::size_type pos(0);
  std::string X=Line.substr(0,MaxLine);
  std::string::size_type posB=X.find_last_of(" ,");
  size_t spc(0);
  while (posB!=std::string::npos && 
	 X.length()>=MaxLine-spc)
    {
      pos+=posB+1;
      if (!isspace(X[posB]))
	posB++;
      const std::string Out=X.substr(0,posB);
      if (!isEmpty(Out))
        {
	  OX<<"c ";
	  if (spc)
	    OX<<std::string(spc,' ');
	  OX<<X.substr(0,posB)<<std::endl;
	}
      spc=8;
      X=Line.substr(pos,MaxLine-spc);
      posB=X.find_last_of(" ,");
    }
  if (!isEmpty(X))
    {
      OX<<"c ";
      if (spc)
	OX<<std::string(spc,' ');
      OX<<X<<std::endl;
    }
  return;
}

}  // NAMESPACE StrFunc