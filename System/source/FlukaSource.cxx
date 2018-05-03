/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   source/FlukaSource.cxx
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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <memory>

#include "Exception.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "support.h"
#include "writeSupport.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "doubleErr.h"
#include "varList.h"
#include "inputSupport.h"
#include "Source.h"
#include "SrcItem.h"
#include "SrcData.h"
#include "surfRegister.h"
#include "ModelSupport.h"
#include "HeadRule.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "FixedOffset.h"
#include "WorkData.h"
#include "World.h"
#include "particleConv.h"

#include "SourceBase.h"
#include "FlukaSource.h"

namespace SDef
{

FlukaSource::FlukaSource(const std::string& keyName) : 
  FixedOffset(keyName,0),SourceBase()
  /*!
    Constructor BUT ALL variable are left unpopulated.
    \param keyName :: main name
  */
{
  for(size_t i=0;i<12;i++)
    sValues[i]=unitTYPE(0,"");
}

FlukaSource::FlukaSource(const FlukaSource& A) : 
  attachSystem::FixedOffset(A),SourceBase(A),
  sourceName(A.sourceName),sValues(A.sValues)
  /*!
    Copy constructor
    \param A :: FlukaSource to copy
  */
{}

FlukaSource&
FlukaSource::operator=(const FlukaSource& A)
  /*!
    Assignment operator
    \param A :: FlukaSource to copy
    \return *this
  */
{
  if (this!=&A)
    {
      attachSystem::FixedOffset::operator=(A);
      SourceBase::operator=(A);
      sourceName=A.sourceName;
      sValues=A.sValues;
    }
  return *this;
}

FlukaSource::~FlukaSource() 
  /*!
    Destructor
  */
{}

FlukaSource*
FlukaSource::clone() const
  /*!
    Clone constructor
    \return copy of this
  */
{
  return new FlukaSource(*this);
}
  
  
  
void
FlukaSource::populate(const ITYPE& inputMap)
  /*!
    Populate Varaibles
    \param inputMap :: Control variables
   */
{
  ELog::RegMethod RegA("FlukaSource","populate");

  attachSystem::FixedOffset::populate(inputMap);
  SourceBase::populate(inputMap);


  const size_t NPts=mainSystem::sizeInput(inputMap,"sourceVar");
  for(size_t index=0;index<NPts && index<12;index++)
    {
      double D;
      const std::string IStr=
	mainSystem::getInput<std::string>(inputMap,"sourceVar",index);
      
      if (StrFunc::convert(IStr,D))
	sValues[index]=unitTYPE(1,IStr);
      else if (!IStr.empty() &&
	       IStr!="-" &&
	       StrFunc::toUpperString(IStr)!="DEF")
	sValues[index]=unitTYPE(-1,IStr);
    }
  
  sourceName=mainSystem::getDefInput<std::string>
			   (inputMap,"sourceName",0,"");
  
  
  return;
}

void
FlukaSource::createUnitVector(const attachSystem::FixedComp& FC,
			      const long int linkIndex)
  /*!
    Create the unit vector
    \param FC :: Fixed Component
    \param linkIndex :: Link index [signed for opposite side]
   */
{
  ELog::RegMethod RegA("FlukaSource","createUnitVector");

  attachSystem::FixedComp::createUnitVector(FC,linkIndex);
  applyOffset();
  return;
}

void
FlukaSource::rotate(const localRotate& LR)
  /*!
    Rotate the source
    \param LR :: Rotation to apply
  */
{
  ELog::RegMethod Rega("FlukaSource","rotate");
  FixedComp::applyRotation(LR);  
  return;
}
  
void
FlukaSource::createSource(SDef::Source& sourceCard) const
  /*!
    Creates a simple beam sampled uniformly in a
    circle
    \param sourceCard :: Source system
  */
{
  ELog::RegMethod RegA("FlukaSource","createSource");

  sourceCard.setComp("vec",Y);
  sourceCard.setComp("axs",Y);
    
  sourceCard.setComp("pos",Origin);
  
  SourceBase::createEnergySource(sourceCard);

  return;
}  

void
FlukaSource::createAll(const attachSystem::FixedComp& FC,
		      const long int linkIndex)
  /*!
    Create all with out using Control variables
    \param FC :: Fixed Point for origin/axis of beam
    \param linkIndex :: link Index				
  */
{
  ELog::RegMethod RegA("FlukaSource","createAll(FC)");
  createUnitVector(FC,linkIndex);
  return;
}

  
void
FlukaSource::createAll(const ITYPE& inputMap,
		      const attachSystem::FixedComp& FC,
		      const long int linkIndex)

  /*!
    Create all the source
    \param Control :: DataBase for variables
    \param FC :: Fixed Point for origin/axis of beam
    \param linkIndex :: link Index				
   */
{
  ELog::RegMethod RegA("FlukaSource","createAll<FC,linkIndex>");
  populate(inputMap);
  createUnitVector(FC,linkIndex);

  return;
}

void
FlukaSource::write(std::ostream& OX) const
  /*!
    Write out as a MCNP source system
    \param OX :: Output stream
  */
{
  ELog::RegMethod RegA("FlukaSource","write");

  Source sourceCard;
  createSource(sourceCard);
  sourceCard.write(OX);
  return;
}

void
FlukaSource::writePHITS(std::ostream& OX) const
  /*!
    Write out as a PHITS source system
    \param OX :: Output stream
  */
{
  ELog::RegMethod RegA("FlukaSource","writePHITS");

  return;
}

void
FlukaSource::writeFLUKA(std::ostream& OX) const
  /*!
    Write out as a FLUKA source system
    \param OX :: Output stream
  */
{
  ELog::RegMethod RegA("FlukaSource","writeFLUKA");

  const particleConv& PC=particleConv::Instance();

  // can be two for an energy range
  if (Energy.size()!=1)
    throw ColErr::SizeError<size_t>
      (Energy.size(),1,"Energy only single point supported");

  std::ostringstream cx;
  // energy : energy divirgence : angle spread [mrad]
  // radius : innerRadius : -1 t o means radius
  cx<<"BEAM "<<-0.001*Energy.front()<<" 0.0 "<<" - "
    <<" "<<" - "<<" - - ";
  cx<<StrFunc::toUpperString(particleType);
  StrFunc::writeFLUKA(cx.str(),OX);
  cx.str("");

  // Y Axis is Z in fluka, X is X
  cx<<"BEAMAXES "<<X<<" "<<Y;
  StrFunc::writeFLUKA(cx.str(),OX);
  cx.str("");
  cx<<"BEAMPOS "<<Origin;
  StrFunc::writeFLUKA(cx.str(),OX);
  cx.str("");
  OX<<"SOURCE"<<std::endl;
  
  return;
}




} // NAMESPACE SDef
