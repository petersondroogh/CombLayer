/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   balder/CrossPipe.cxx
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
#include <array>

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
#include "Quadratic.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Rules.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "Simulation.h"
#include "ModelSupport.h"
#include "MaterialSupport.h"
#include "generateSurf.h"
#include "LinkUnit.h"  
#include "FixedComp.h"
#include "FixedOffset.h"
#include "ContainedComp.h"
#include "BaseMap.h"
#include "CellMap.h"
#include "SurfMap.h"
#include "FrontBackCut.h"
#include "SurfMap.h"
#include "SurInter.h"

#include "CrossPipe.h"

namespace xraySystem
{

CrossPipe::CrossPipe(const std::string& Key) : 
  attachSystem::FixedOffset(Key,6),
  attachSystem::ContainedComp(),attachSystem::CellMap(),
  attachSystem::SurfMap(),attachSystem::FrontBackCut(),
  vacIndex(ModelSupport::objectRegister::Instance().cell(Key)),
  cellIndex(vacIndex+1)
  /*!
    Constructor BUT ALL variable are left unpopulated.
    \param Key :: KeyName
  */
{}


CrossPipe::~CrossPipe() 
  /*!
    Destructor
  */
{}

void
CrossPipe::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: DataBase of variables
  */
{
  ELog::RegMethod RegA("CrossPipe","populate");
  
  FixedOffset::populate(Control);

  // Void + Fe special:
  horrRadius=Control.EvalVar<double>(keyName+"HorrRadius");
  vertRadius=Control.EvalVar<double>(keyName+"VertRadius");
  height=Control.EvalVar<double>(keyName+"Height");
  depth=Control.EvalVar<double>(keyName+"Depth");
  
  length=Control.EvalVar<double>(keyName+"Length");

  feThick=Control.EvalVar<double>(keyName+"FeThick");
  topPlate=Control.EvalVar<double>(keyName+"TopPlate");
  basePlate=Control.EvalVar<double>(keyName+"BasePlate");
  
  flangeRadius=Control.EvalDefVar<double>(keyName+"FlangeRadius",-1.0);
  flangeLength=Control.EvalDefVar<double>(keyName+"FlangeLength",0.0);
  
  voidMat=ModelSupport::EvalDefMat<int>(Control,keyName+"VoidMat",0);
  feMat=ModelSupport::EvalMat<int>(Control,keyName+"FeMat");

  return;
}

void
CrossPipe::createUnitVector(const attachSystem::FixedComp& FC,
                             const long int sideIndex)
  /*!
    Create the unit vectors
    \param FC :: Fixed component to link to
    \param sideIndex :: Link point and direction [0 for origin]
  */
{
  ELog::RegMethod RegA("CrossPipe","createUnitVector");

  FixedComp::createUnitVector(FC,sideIndex);
  yStep+=length/2.0;
  applyOffset();

  return;
}

  
void
CrossPipe::getShiftedSurf(const HeadRule& HR,
			   const int index,
			   const int dFlag,
			   const double length)
  /*!
    Support function to calculate the shifted surface
    \param HR :: HeadRule to extract plane surf
    \param index :: offset index
    \param dFlag :: direction flag
    \param length :: length to shift by
  */
{
  ELog::RegMethod RegA("CrossPipe","getShiftedSurf");
  
  std::set<int> FS=HR.getSurfSet();
  for(const int& SN : FS)
    {
      const Geometry::Surface* SPtr=SMap.realSurfPtr(SN);

      const Geometry::Plane* PPtr=
	dynamic_cast<const Geometry::Plane*>(SPtr);
      if (PPtr)
	{
	  if (SN*dFlag>0)
	    ModelSupport::buildShiftedPlane
	      (SMap,vacIndex+index,PPtr,dFlag*length);
	  else
	    ModelSupport::buildShiftedPlaneReversed
	      (SMap,vacIndex+index,PPtr,dFlag*length);
	  
	  return;
	}
      const Geometry::Cylinder* CPtr=
	dynamic_cast<const Geometry::Cylinder*>(SPtr);
      // Cylinder case:
      if (CPtr)
	{
	  if (SN>0)
	    ModelSupport::buildCylinder
	      (SMap,vacIndex+index,CPtr->getCentre()+Y*length,
	       CPtr->getNormal(),CPtr->getRadius());
	  else
	    ModelSupport::buildCylinder
	      (SMap,vacIndex+index,CPtr->getCentre()-Y*flangeLength,
	       CPtr->getNormal(),CPtr->getRadius());
	  return;
	}
    }
  throw ColErr::EmptyValue<int>("HeadRule contains no planes/cylinder");
} 

void
CrossPipe::createSurfaces()
  /*!
    Create the surfaces
  */
{
  ELog::RegMethod RegA("CrossPipe","createSurfaces");
  
  // middle of the flange
  //  const double midFlange((length-flangeLength)/2.0);
  
  // Inner void
  if (frontActive())
    // create surface 101:
    getShiftedSurf(getFrontRule(),101,1,flangeLength);
  else
    {
      ModelSupport::buildPlane(SMap,vacIndex+1,Origin-Y*(length/2.0),Y);
      ModelSupport::buildPlane(SMap,vacIndex+101,Origin-Y*(length/2.0-flangeLength),Y);
      FrontBackCut::setFront(SMap.realSurf(vacIndex+1));
    }
  
  // Inner void
  if (backActive())
    // create surface 102:
    getShiftedSurf(getBackRule(),102,-1,flangeLength);
  else
    {
      ModelSupport::buildPlane(SMap,vacIndex+2,Origin+Y*(length/2.0),Y);
      ModelSupport::buildPlane(SMap,vacIndex+102,
			       Origin+Y*(length/2.0-flangeLength),Y);
      FrontBackCut::setBack(-SMap.realSurf(vacIndex+2));
    }

  
  // MAIN SURFACES axial:
  ModelSupport::buildCylinder(SMap,vacIndex+7,Origin,Y,horrRadius);
  ModelSupport::buildCylinder(SMap,vacIndex+17,Origin,Y,horrRadius+feThick);

  // FLANGE SURFACES:
  if (flangeRadius>Geometry::zeroTol && flangeLength>Geometry::zeroTol)
    ModelSupport::buildCylinder(SMap,vacIndex+107,Origin,Y,horrRadius+flangeRadius);

  // Secondary SURFACES:
  ModelSupport::buildCylinder(SMap,vacIndex+207,Origin,Z,vertRadius);
  ModelSupport::buildCylinder(SMap,vacIndex+217,Origin,Z,vertRadius+feThick);

  // top plates
  ModelSupport::buildPlane(SMap,vacIndex+206,Origin+Z*height,Z);
  ModelSupport::buildPlane(SMap,vacIndex+216,Origin+Z*(height+topPlate),Z);

  // base plates
  ModelSupport::buildPlane(SMap,vacIndex+205,Origin-Z*depth,Z);
  ModelSupport::buildPlane(SMap,vacIndex+215,Origin-Z*(depth+basePlate),Z);

  return;
}

void
CrossPipe::createObjects(Simulation& System)
  /*!
    Adds the vacuum box
    \param System :: Simulation to create objects in
   */
{
  ELog::RegMethod RegA("CrossPipe","createObjects");

  std::string Out;
  
  const std::string frontStr=frontRule();
  const std::string backStr=backRule();
  
  // Void [want a single void]
  Out=ModelSupport::getComposite(SMap,vacIndex," (-7 : -207) 205 -206  ");
  System.addCell(MonteCarlo::Qhull(cellIndex++,voidMat,0.0,Out+frontStr+backStr));
  addCell("Void",cellIndex-1);

  // Main steel pipe
  Out=ModelSupport::getComposite(SMap,vacIndex," -17 7 207 ");
  System.addCell(MonteCarlo::Qhull(cellIndex++,feMat,0.0,Out+frontStr+backStr));
  addCell("Pipe",cellIndex-1);

  // Vert steel pipe
  Out=ModelSupport::getComposite(SMap,vacIndex," 17 -217 207 205 -206 ");
  System.addCell(MonteCarlo::Qhull(cellIndex++,feMat,0.0,Out));
  addCell("Pipe",cellIndex-1);

  // BasePlate
  Out=ModelSupport::getComposite(SMap,vacIndex," -217 215 -205 ");
  System.addCell(MonteCarlo::Qhull(cellIndex++,feMat,0.0,Out));
  addCell("BasePlate",cellIndex-1);

  // BasePlate
  Out=ModelSupport::getComposite(SMap,vacIndex," -217 206 -216 ");
  System.addCell(MonteCarlo::Qhull(cellIndex++,feMat,0.0,Out));
  addCell("TopPlate",cellIndex-1);
  

  if (flangeRadius>Geometry::zeroTol && flangeLength>Geometry::zeroTol)
    {
      Out=ModelSupport::getComposite(SMap,vacIndex," 17 -107 -101 ");
      System.addCell(MonteCarlo::Qhull(cellIndex++,feMat,0.0,Out+frontStr));
      addCell("Flange",cellIndex-1);

      Out=ModelSupport::getComposite(SMap,vacIndex," 17 -107 102 ");
      System.addCell(MonteCarlo::Qhull(cellIndex++,feMat,0.0,Out+backStr));
      addCell("Flange",cellIndex-1);

      Out=ModelSupport::getComposite(SMap,vacIndex," 17 -107 101 -102 217 ");
      System.addCell(MonteCarlo::Qhull(cellIndex++,0,0.0,Out));
      addCell("FlangeVoid",cellIndex-1);
      Out=ModelSupport::getComposite(SMap,vacIndex," (-107 : -217) 215 -216  ");
    }
  else
    Out=ModelSupport::getComposite(SMap,vacIndex," (-17 : -217) 215 -216  ");


  addOuterSurf(Out+frontStr+backStr);

  return;
}
  
void
CrossPipe::createLinks()
  /*!
    Determines the link point on the outgoing plane.
    It must follow the beamline, but exit at the plane
  */
{
  ELog::RegMethod RegA("CrossPipe","createLinks");

  //stufff for intersection


  FrontBackCut::createLinks(*this,Origin,Y);  //front and back

  FixedComp::setConnect(2,Origin-X*horrRadius,-X);
  FixedComp::setConnect(3,Origin+X*horrRadius,X);
  FixedComp::setConnect(4,Origin-Z*depth,-Z);
  FixedComp::setConnect(5,Origin+Z*height,Z);
  
  FixedComp::setLinkSurf(2,SMap.realSurf(vacIndex+7));
  FixedComp::setLinkSurf(3,SMap.realSurf(vacIndex+7));
  FixedComp::setLinkSurf(4,SMap.realSurf(vacIndex+5));
  FixedComp::setLinkSurf(5,SMap.realSurf(vacIndex+6));

  /*
  FixedComp::setConnect(7,Origin-Z*(radius+feThick),-Z);
  FixedComp::setConnect(8,Origin+Z*(radius+feThick),Z);
  FixedComp::setLinkSurf(7,SMap.realSurf(vacIndex+17));
  FixedComp::setLinkSurf(8,SMap.realSurf(vacIndex+17));
      
  
  // MID Point: [NO SURF]
  const Geometry::Vec3D midPt=
    (getLinkPt(1)+getLinkPt(2))/2.0;
  FixedComp::setConnect(6,midPt,Y);

  if (flangeRadius>0.0)
    {
      FixedComp::setConnect(9,Origin-Z*horrRadius+flangeRadius,-Z);
      FixedComp::setConnect(10,Origin+Z*horrRadius+flangeRadius,Z);
      
      FixedComp::setLinkSurf(9,SMap.realSurf(vacIndex+107));
      FixedComp::setLinkSurf(10,SMap.realSurf(vacIndex+107));
    }
  else
    {
      FixedComp::setConnect(9,Origin-Z*(flangeHeight/2.0),-Z);
      FixedComp::setConnect(10,Origin+Z*(flangeHeight/2.0),Z);
      FixedComp::setLinkSurf(9,-SMap.realSurf(vacIndex+105));
      FixedComp::setLinkSurf(10,SMap.realSurf(vacIndex+106));
    }
  */
  return;
}
  
  
void
CrossPipe::createAll(Simulation& System,
		      const attachSystem::FixedComp& FC,
		      const long int FIndex)
 /*!
    Generic function to create everything
    \param System :: Simulation item
    \param FC :: FixedComp
    \param FIndex :: Fixed Index
  */
{
  ELog::RegMethod RegA("CrossPipe","createAll(FC)");

  populate(System.getDataBase());
  createUnitVector(FC,FIndex);
  createSurfaces();    
  createObjects(System);
  createLinks();
  insertObjects(System);   
  
  return;
}
  
}  // NAMESPACE xraySystem
