/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   essBuild/RefPlug.cxx
*
 * Copyright (c) 2004-2013 by Stuart Ansell
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
#include <boost/shared_ptr.hpp>

#include "Exception.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Quaternion.h"
#include "Surface.h"
#include "surfIndex.h"
#include "Quadratic.h"
#include "Rules.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "Simulation.h"
#include "ModelSupport.h"
#include "generateSurf.h"
#include "support.h"
#include "stringCombine.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "ContainedComp.h"
#include "LayerComp.h"
#include "RefPlug.h"

namespace snsSystem
{

RefPlug::RefPlug(const std::string& Key) :
  attachSystem::ContainedComp(),attachSystem::LayerComp(),
  attachSystem::FixedComp(Key,3),
  refIndex(ModelSupport::objectRegister::Instance().cell(Key)),
  cellIndex(refIndex+1)
  /*!
    Constructor
    \param Key :: Name of construction key
  */
{}

RefPlug::RefPlug(const RefPlug& A) : 
  attachSystem::ContainedComp(A),attachSystem::LayerComp(A),
  attachSystem::FixedComp(A), 
  refIndex(A.refIndex),cellIndex(A.cellIndex),xStep(A.xStep),
  yStep(A.yStep),zStep(A.zStep),xyAngle(A.xyAngle),
  zAngle(A.zAngle),nLayers(A.nLayers),radius(A.radius),
  height(A.height),mat(A.mat),temp(A.temp)
  /*!
    Copy constructor
    \param A :: RefPlug to copy
  */
{}

RefPlug&
RefPlug::operator=(const RefPlug& A)
  /*!
    Assignment operator
    \param A :: RefPlug to copy
    \return *this
  */
{
  if (this!=&A)
    {
      attachSystem::ContainedComp::operator=(A);
      attachSystem::LayerComp::operator=(A);
      attachSystem::FixedComp::operator=(A);
      cellIndex=A.cellIndex;
      xStep=A.xStep;
      yStep=A.yStep;
      zStep=A.zStep;
      xyAngle=A.xyAngle;
      zAngle=A.zAngle;
      nLayers=A.nLayers;
      radius=A.radius;
      height=A.height;
      mat=A.mat;
      temp=A.temp;
    }
  return *this;
}

RefPlug::~RefPlug()
  /*!
    Destructor
  */
{}

RefPlug*
RefPlug::clone() const
  /*!
    Clone copy constructor
    \return copy of this
  */
{
  return new RefPlug(*this);
}

void
RefPlug::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable table to use
  */
{
  ELog::RegMethod RegA("RefPlug","populate");

    // Master values
  xStep=Control.EvalVar<double>(keyName+"XStep");
  yStep=Control.EvalVar<double>(keyName+"YStep");
  zStep=Control.EvalVar<double>(keyName+"ZStep");
  xyAngle=Control.EvalVar<double>(keyName+"XYangle");
  zAngle=Control.EvalVar<double>(keyName+"Zangle");

  double R,H,T;
  int M;
  nLayers=Control.EvalVar<size_t>(keyName+"NLayers");   
  for(size_t i=0;i<=nLayers;i++)
    {
      if (i)
	{
	  H+=2.0*Control.EvalVar<double>
	    (StrFunc::makeString(keyName+"HGap",i));   
	  R+=Control.EvalVar<double>
	    (StrFunc::makeString(keyName+"RadGap",i));   
	  M=Control.EvalVar<int>
	    (StrFunc::makeString(keyName+"Material",i));   
	  T=(!M) ? 0.0 : 
	    Control.EvalVar<double>(StrFunc::makeString(keyName+"Temp",i)); 
	}
      else
	{
	  H=Control.EvalVar<double>(keyName+"Height");   
	  R=Control.EvalVar<double>(keyName+"Radius");   
	  M=Control.EvalVar<int>(keyName+"Mat");   
	  T=Control.EvalVar<int>(keyName+"Temp");   
	}
      radius.push_back(R);
      height.push_back(H);
      mat.push_back(M);
      temp.push_back(T);
    }

  return;
}

void
RefPlug::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: Fixed Component
  */
{
  ELog::RegMethod RegA("RefPlug","createUnitVector");
  attachSystem::FixedComp::createUnitVector(FC);
  applyShift(xStep,yStep,zStep);
  applyAngleRotate(xyAngle,zAngle);

  return;
}

void
RefPlug::createSurfaces()
  /*!
    Create planes for the silicon and Polyethene layers
  */
{
  ELog::RegMethod RegA("RefPlug","createSurfaces");

  // Divide plane
  ModelSupport::buildPlane(SMap,refIndex+1,Origin,X);  
  ModelSupport::buildPlane(SMap,refIndex+2,Origin,Y);  

  int SI(refIndex);
  for(size_t i=0;i<=nLayers;i++)
    {
      ModelSupport::buildCylinder(SMap,SI+7,Origin,Z,radius[i]);  
      ModelSupport::buildPlane(SMap,SI+5,Origin-Z*height[i]/2.0,Z);  
      ModelSupport::buildPlane(SMap,SI+6,Origin+Z*height[i]/2.0,Z);  
      SI+=10;
    }

  // Across sides
  return; 
}

void
RefPlug::createObjects(Simulation& System)
  /*!
    Create the vaned moderator
    \param System :: Simulation to add results
  */
{
  ELog::RegMethod RegA("RefPlug","createObjects");

  std::string Out;
  int SI(refIndex);
  for(size_t i=0;i<=nLayers;i++)
    {
      Out=ModelSupport::getComposite(SMap,SI," -7 5 -6 ");

      if (i==nLayers) addOuterSurf(Out);
      if (i)
	Out+=ModelSupport::getComposite(SMap,SI-10," (7:-5:6) ");
      System.addCell(MonteCarlo::Qhull(cellIndex++,mat[i],temp[i],Out));
      SI+=10;
    }
  return; 
}

void
RefPlug::createLinks()
  /*!
    Creates a full attachment set
  */
{  
  ELog::RegMethod RegA("RefPlug","createLinks");

  if (nLayers)
    {
      const int SI(refIndex+static_cast<int>(nLayers)*10);

      FixedComp::setConnect(0,Origin-Y*radius[nLayers],-Y);
      FixedComp::setLinkSurf(0,SMap.realSurf(SI+7));
      FixedComp::addLinkSurf(0,-SMap.realSurf(refIndex+2));

      FixedComp::setConnect(1,Origin+Y*radius[nLayers],Y);
      FixedComp::setLinkSurf(1,SMap.realSurf(SI+7));
      FixedComp::addLinkSurf(1,SMap.realSurf(refIndex+2));

      FixedComp::setConnect(2,Origin-X*radius[nLayers],-X);
      FixedComp::setLinkSurf(2,SMap.realSurf(SI+7));
      FixedComp::addLinkSurf(2,-SMap.realSurf(refIndex+1));

      FixedComp::setConnect(3,Origin+X*radius[nLayers],X);
      FixedComp::setLinkSurf(3,SMap.realSurf(SI+7));
      FixedComp::addLinkSurf(3,SMap.realSurf(refIndex+1));
      
      FixedComp::setConnect(4,Origin-Z*(height[nLayers]/2.0),-Z);
      FixedComp::setLinkSurf(4,-SMap.realSurf(SI+5));

      FixedComp::setConnect(5,Origin+Z*(height[nLayers]/2.0),Z);
      FixedComp::setLinkSurf(5,SMap.realSurf(SI+6));
    }
  else 
    ELog::EM<<"NO Layers in RefPlug"<<ELog::endErr;
  return;
}



Geometry::Vec3D
RefPlug::getSurfacePoint(const size_t layerIndex,
			      const size_t sideIndex) const
  /*!
    Given a side and a layer calculate the link point
    \param sideIndex :: Side [0-5]
    \param layerIndex :: layer, 0 is inner moderator [0-6]
    \return Surface point
  */
{
  ELog::RegMethod RegA("RefPlug","getSurfacePoint");

  if (sideIndex>5) 
    throw ColErr::IndexError<size_t>(sideIndex,5,"sideIndex ");
  if (layerIndex>=nLayers) 
    throw ColErr::IndexError<size_t>(layerIndex,nLayers,"layer");

  // Modification map:
  switch(sideIndex)
    {
    case 0:
      return Origin-Y*radius[layerIndex];
    case 1:
      return Origin+Y*radius[layerIndex];
    case 2:
      return Origin-X*radius[layerIndex];
    case 3:
      return Origin+X*radius[layerIndex];
    case 4:
      return Origin-Z*(height[layerIndex]/2.0);
    case 5:
      return Origin+Z*(height[layerIndex]/2.0);
    }
  throw ColErr::IndexError<size_t>(sideIndex,5,"sideIndex ");
}

std::string
RefPlug::getLayerString(const size_t layerIndex,
			     const size_t sideIndex) const
  /*!
    Given a side and a layer calculate the link surf
    \param sideIndex :: Side [0-5]
    \param layerIndex :: layer, 0 is inner moderator [0-4]
    \return Surface string
  */
{
  ELog::RegMethod RegA("RefPlug","getLayerString");

  if (layerIndex>=nLayers) 
    throw ColErr::IndexError<size_t>(layerIndex,nLayers,"layer");

  const int SI(refIndex+static_cast<int>(layerIndex)*10);
  std::ostringstream cx;
  switch(sideIndex)
    {
    case 0:
      return ModelSupport::getComposite(SMap,SI,refIndex," 7 -2M ");
    case 1:
      return ModelSupport::getComposite(SMap,SI,refIndex," 7 2M ");
    case 2:
      return ModelSupport::getComposite(SMap,SI,refIndex," 7 -1M ");
    case 3:
      return ModelSupport::getComposite(SMap,SI,refIndex," 7 1M ");
    case 4:
      cx<<" "<<-SMap.realSurf(SI+5)<<" ";
      return cx.str();
    case 5:
      cx<<" "<<SMap.realSurf(SI+6)<<" ";
      return cx.str();
    }
  throw ColErr::IndexError<size_t>(sideIndex,5,"sideIndex ");
}

int
RefPlug::getLayerSurf(const size_t layerIndex,
			 const size_t sideIndex) const
  /*!
    Given a side and a layer calculate the link surf
    \param sideIndex :: Side [0-5]
    \param layerIndex :: layer, 0 is inner moderator [0-4]
    \return Surface string
  */
{
  ELog::RegMethod RegA("RefPlug","getLayerSurf");

  if (layerIndex>=nLayers) 
    throw ColErr::IndexError<size_t>(layerIndex,nLayers,"layer");
  
  const int SI(refIndex+static_cast<int>(layerIndex)*10);
  switch(sideIndex)
    {
    case 0:
      return SMap.realSurf(SI+7);
    case 1:
      return SMap.realSurf(SI+7);
    case 2:
      return SMap.realSurf(SI+7);
    case 3:
      return SMap.realSurf(SI+7);
    case 4:
      return -SMap.realSurf(SI+5);
    case 5:
      return SMap.realSurf(SI+6);
    }
  throw ColErr::IndexError<size_t>(sideIndex,5,"sideIndex ");
}

int
RefPlug::getCommonSurf(const size_t sideIndex) const
  /*!
    Given a side calculate the boundary surface
    \param sideIndex :: Side [0-5]
    \return Common dividing surface [outward pointing]
  */
{
  ELog::RegMethod RegA("RefPlug","getCommonSurf");

  switch(sideIndex)
    {
    case 0:
      return -SMap.realSurf(refIndex+2);
    case 1:
      return SMap.realSurf(refIndex+2);
    case 2:
      return -SMap.realSurf(refIndex+1);
    case 3:
      return SMap.realSurf(refIndex+1);
    case 4:
    case 5:
      return 0;
    }
  throw ColErr::IndexError<size_t>(sideIndex,5,"sideIndex ");
}

void
RefPlug::createAll(Simulation& System,
		     const attachSystem::FixedComp& FC)
  /*!
    Extrenal build everything
    \param System :: Simulation
    \param FC :: FixedComponent for origin
   */
{
  ELog::RegMethod RegA("RefPlug","createAll");
  populate(System.getDataBase());

  createUnitVector(FC);
  createSurfaces();
  createObjects(System);
  createLinks();
  insertObjects(System);       

  return;
}

}  // NAMESPACE instrumentSystem
