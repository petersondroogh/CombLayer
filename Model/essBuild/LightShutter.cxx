/*********************************************************************
  CombLayer : MCNP(X) Input builder

 * File:   essBuild/LightShutter.cxx
 *
 * Copyright (c) 2004-2017 by Stuart Ansell
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
#include "stringCombine.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Quaternion.h"
#include "Surface.h"
#include "surfIndex.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "surfEqual.h"
#include "Quadratic.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Line.h"
#include "Rules.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "inputParam.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "Simulation.h"
#include "ReadFunctions.h"
#include "ModelSupport.h"
#include "MaterialSupport.h"
#include "generateSurf.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "FixedOffset.h"
#include "ContainedComp.h"
#include "BaseMap.h"
#include "CellMap.h"
#include "surfDBase.h"
#include "surfDIter.h"
#include "surfDivide.h"
#include "SurInter.h"
#include "mergeTemplate.h"

#include "LightShutter.h"

namespace essSystem
{

LightShutter::LightShutter(const std::string& Key)  :
  attachSystem::ContainedComp(),
  attachSystem::FixedOffset(Key,8),
  attachSystem::CellMap(),
  lightIndex(ModelSupport::objectRegister::Instance().cell(Key)),
  cellIndex(lightIndex+1)
  /*!
    Constructor BUT ALL variable are left unpopulated.
    \param Key :: Name for item in search
  */
{}

LightShutter::LightShutter(const LightShutter& A) :
  attachSystem::ContainedComp(A),attachSystem::FixedOffset(A),
  attachSystem::CellMap(A),
  lightIndex(A.lightIndex),cellIndex(A.cellIndex),
  active(A.active),
  length(A.length),width(A.width),height(A.height),
  nLayers(A.nLayers),
  wallThick(A.wallThick),mainMat(A.mainMat),wallMat(A.wallMat)
  /*!
    Copy constructor
    \param A :: LightShutter to copy
  */
{}

LightShutter&
LightShutter::operator=(const LightShutter& A)
  /*!
    Assignment operator
    \param A :: LightShutter to copy
    \return *this
  */
{
  if (this!=&A)
    {
      attachSystem::ContainedComp::operator=(A);
      attachSystem::FixedOffset::operator=(A);
      attachSystem::CellMap::operator=(A);
      cellIndex=A.cellIndex;
      active=A.active;
      length=A.length;
      width=A.width;
      height=A.height;
      nLayers=A.nLayers;
      wallThick=A.wallThick;
      mainMat=A.mainMat;
      wallMat=A.wallMat;
    }
  return *this;
}

LightShutter*
LightShutter::clone() const
  /*!
    virtual copy constructor
    \return new LightShutter(*this)
  */
{
  return new LightShutter(*this);
}

LightShutter::~LightShutter()
  /*!
    Destructor
  */
{}

  void
  LightShutter::layerProcess(Simulation& System, const std::string& cellName,
			     const size_t& lpS, const size_t& lsS,
			     const int& nLayers, const int& mat)
  /*!
    Processes the splitting of the surfaces into a multilayer system
    \param System :: Simulation to work on
    \param cellName :: Name of the cell to split
    \param lpS :: link pont of primary surface
    \param lsS :: link point of secondary surface
    \param nLayers :: number of layers to divide to
    \param mat :: material
  */
  {
    ELog::RegMethod RegA("LightShutter","layerProcess");
    if (nLayers>1)
      {
	const int pS = getLinkSurf(lpS);
	const int sS = getLinkSurf(lsS);

	const attachSystem::CellMap* CM = dynamic_cast<const attachSystem::CellMap*>(this);
	MonteCarlo::Object* wallObj(0);
	int wallCell(0);

	if (CM)
	  {
	    wallCell=CM->getCell(cellName);
	    wallObj=System.findQhull(wallCell);
	  }

	if (!wallObj)
	  throw ColErr::InContainerError<int>(wallCell,
					      "Cell '" + cellName + "' not found");

	double baseFrac = 1.0/nLayers;
	ModelSupport::surfDivide DA;
	for(int i=1;i<nLayers;i++)
	  {
	    DA.addFrac(baseFrac);
	    DA.addMaterial(mat);
	    baseFrac += 1.0/nLayers;
	  }
	DA.addMaterial(mat);

	DA.setCellN(wallCell);
	DA.setOutNum(cellIndex, lightIndex+10000);

	ModelSupport::mergeTemplate<Geometry::Plane,
				    Geometry::Plane> surroundRule;

	surroundRule.setSurfPair(SMap.realSurf(pS),
				 SMap.realSurf(sS));

	std::string OutA = getLinkString(lpS);
	std::string OutB = getLinkComplement(lsS);

	surroundRule.setInnerRule(OutA);
	surroundRule.setOuterRule(OutB);

	DA.addRule(&surroundRule);
	DA.activeDivideTemplate(System);

	cellIndex=DA.getCellNum();
      }
  }

  
void
LightShutter::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable data base
  */
{
  ELog::RegMethod RegA("LightShutter","populate");

  FixedOffset::populate(Control);

  active=Control.EvalDefVar<int>(keyName+"Active", 0);
  if (!active) return;

  height=Control.EvalVar<double>(keyName+"Height");
  nLayers=Control.EvalVar<int>(keyName+"NLayers");
  width=Control.EvalVar<double>(keyName+"Width");
  length=Control.EvalVar<double>(keyName+"Length");
  wallThick=Control.EvalVar<double>(keyName+"WallThick");

  mainMat=ModelSupport::EvalMat<int>(Control,keyName+"MainMat");
  wallMat=ModelSupport::EvalMat<int>(Control,keyName+"WallMat");

  return;
}

void
LightShutter::createUnitVector(const attachSystem::FixedComp& FC,
			       const long int sideIndex)

  /*!
    Create the unit vectors
    \param FC :: object for origin/radius
    \param sideIndex :: Side for monolith
  */
{
  ELog::RegMethod RegA("LightShutter","createUnitVector");
  FixedComp::createUnitVector(FC,sideIndex);
  applyOffset();
  return;
}

void
LightShutter::createSurfaces()
  /*!
    Create All the surfaces
  */
{
  ELog::RegMethod RegA("LightShutter","createSurface");


  ModelSupport::buildPlane(SMap,lightIndex+1,Origin,Y);
  ModelSupport::buildPlane(SMap,lightIndex+2,Origin+Y*length,Y);
  ModelSupport::buildPlane(SMap,lightIndex+3,Origin-X*(width/2.0),X);
  ModelSupport::buildPlane(SMap,lightIndex+4,Origin+X*(width/2.0),X);
  ModelSupport::buildPlane(SMap,lightIndex+5,Origin-Z*(height/2.0),Z);
  ModelSupport::buildPlane(SMap,lightIndex+6,Origin+Z*(height/2.0),Z);

  ModelSupport::buildPlane(SMap,lightIndex+12,
  			   Origin+Y*(length+wallThick),Y);
  ModelSupport::buildPlane(SMap,lightIndex+15,
			   Origin-Z*(height/2.0+wallThick),Z);

  return;
}


void
LightShutter::createObjects(Simulation& System)
  /*!
    Adds the all the components
    \param System :: Simulation to create objects in
  */
{
  ELog::RegMethod RegA("LightShutter","createObjects");

  if (!active) return;

  std::string Out;
  // Tugnsten middle
  Out=ModelSupport::getComposite(SMap,lightIndex," 1 -2 3 -4 5 -6 ");
  System.addCell(MonteCarlo::Qhull(cellIndex++,mainMat,0.0,Out));
  setCell("main",cellIndex-1);

  layerProcess(System, "main", 6,7, nLayers, mainMat);

  // HDPE boron sheets
  Out=ModelSupport::getComposite(SMap,lightIndex," 2 -12 3 -4 5 -6 ");
  System.addCell(MonteCarlo::Qhull(cellIndex++,wallMat,0.0,Out));
  setCell("wallBack",cellIndex-1);

  Out=ModelSupport::getComposite(SMap,lightIndex," 1 -12 3 -4 15 -5 ");
  System.addCell(MonteCarlo::Qhull(cellIndex++,wallMat,0.0,Out));
  setCell("wallBottom",cellIndex-1);

  Out=ModelSupport::getComposite(SMap,lightIndex," 1 -12 3 -4 15 -6");
  addOuterSurf(Out);

  return;
}


void
LightShutter::createLinks()
  /*!
    Create all the linkes [OutGoing]
  */
{
  ELog::RegMethod RegA("LightShutter","createLinks");

  FixedComp::setConnect(0,Origin,-Y);
  FixedComp::setConnect(1,Origin+Y*(length+wallThick),Y);
  FixedComp::setConnect(2,Origin-X*(width/2.0),-X);
  FixedComp::setConnect(3,Origin+X*(width/2.0),X);
  FixedComp::setConnect(4,Origin-Z*(wallThick+height/2.0),-Z);
  FixedComp::setConnect(5,Origin+Z*(height/2.0),Z);
  FixedComp::setConnect(6,Origin,Y); // for layerProcess
  FixedComp::setConnect(7,Origin+Y*length,Y); // for layerProcess

  FixedComp::setLinkSurf(0,-SMap.realSurf(lightIndex+1));
  FixedComp::setLinkSurf(1,SMap.realSurf(lightIndex+12));
  FixedComp::setLinkSurf(2,-SMap.realSurf(lightIndex+3));
  FixedComp::setLinkSurf(3,SMap.realSurf(lightIndex+4));
  FixedComp::setLinkSurf(4,-SMap.realSurf(lightIndex+15));
  FixedComp::setLinkSurf(5,SMap.realSurf(lightIndex+6));
  FixedComp::setLinkSurf(6,SMap.realSurf(lightIndex+1)); // for layerProcess
  FixedComp::setLinkSurf(7,SMap.realSurf(lightIndex+2)); // for layerProcess

  return;
}

void
LightShutter::createAll(Simulation& System,
		   const attachSystem::FixedComp& FC,
			const long int sideIndex)
  /*!
    Generic function to create everything
    \param System :: Simulation item
    \param FC :: Central origin
    \param sideIndex :: Side of monolith
  */
{
  ELog::RegMethod RegA("LightShutter","createAll");

  populate(System.getDataBase());
  createUnitVector(FC,sideIndex);
  createSurfaces();
  createLinks();
  createObjects(System);
  insertObjects(System);

  return;
}

}  // NAMESPACE essSystem
