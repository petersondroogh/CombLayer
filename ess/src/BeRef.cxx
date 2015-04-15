/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   essBuild/BeRef.cxx
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
#include "MaterialSupport.h"
#include "generateSurf.h"
#include "support.h"
#include "stringCombine.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "ContainedComp.h"
#include "BeRef.h"


namespace essSystem
{

BeRef::BeRef(const std::string& Key) :
  attachSystem::ContainedComp(),attachSystem::FixedComp(Key,6),
  refIndex(ModelSupport::objectRegister::Instance().cell(Key)),
  cellIndex(refIndex+1)
  /*!
    Constructor
    \param Key :: Name of construction key
  */
{}


BeRef::~BeRef()
  /*!
    Destructor
   */
{}

void
BeRef::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable table to use
  */
{
  ELog::RegMethod RegA("BeRef","populate");


    // Master values
  xStep=Control.EvalVar<double>(keyName+"XStep");
  yStep=Control.EvalVar<double>(keyName+"YStep");
  zStep=Control.EvalVar<double>(keyName+"ZStep");
  xyAngle=Control.EvalVar<double>(keyName+"XYangle");
  zAngle=Control.EvalVar<double>(keyName+"Zangle");
  radius=Control.EvalVar<double>(keyName+"Radius");   
  height=Control.EvalVar<double>(keyName+"Height");   
  wallThick=Control.EvalVar<double>(keyName+"WallThick");   
  voidThick=Control.EvalVar<double>(keyName+"VoidThick");   

  refMat=ModelSupport::EvalMat<int>(Control,keyName+"RefMat");
  //  refMatTop=ModelSupport::EvalMat<int>(Control,keyName+"RefMatTop");
  wallMat=ModelSupport::EvalMat<int>(Control,keyName+"WallMat");   

  cellCooling = Control.EvalDefVar<int>(keyName+"CellCooling", 0);
  nRings = Control.EvalDefVar<int>(keyName+"NRings", 4);
  nLayers = Control.EvalDefVar<int>(keyName+"NLayers", 4.0);
  cellWallThick = Control.EvalDefVar<double>(keyName+"CellWallThick", 0.3);
  cellWallMat = ModelSupport::EvalDefMat<int>(Control,keyName+"CellWallMat", 13000);
  cellCoolThick = Control.EvalDefVar<double>(keyName+"CellCoolThick", 0.5);
  cellCoolMat = ModelSupport::EvalDefMat<int>(Control,keyName+"CellCoolMat", 1011);   

  width  = Control.EvalDefVar<double>(keyName+"Width",  -1);
  length = Control.EvalDefVar<double>(keyName+"Length", -1);
  refMat1=ModelSupport::EvalDefMat<int>(Control,keyName+"RefMat1", 0);

  return;
}

void
BeRef::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: Fixed Component
  */
{
  ELog::RegMethod RegA("BeRef","createUnitVector");
  attachSystem::FixedComp::createUnitVector(FC);
  applyShift(xStep,yStep,zStep);
  applyAngleRotate(xyAngle,zAngle);
  
  return;
}

void
BeRef::createSurfaces()
  /*!
    Create Surfaces for the Be
  */
{
  ELog::RegMethod RegA("BeRef","createSurfaces");
      
  ModelSupport::buildCylinder(SMap,refIndex+7,Origin,Z,radius);   // SMap - local surface map [2:1300]
  ModelSupport::buildCylinder(SMap,refIndex+17,Origin,Z,radius+wallThick);  
  ModelSupport::buildCylinder(SMap,refIndex+18,Origin,Z,radius+wallThick+voidThick);  
  
  ModelSupport::buildPlane(SMap,refIndex+6,Origin+Z*(height),Z);  
  ModelSupport::buildPlane(SMap,refIndex+15, Origin+Z*(wallThick),Z); // upper/lower plane of bottom Al container (other plane is the target upper/lower surface)
  ModelSupport::buildPlane(SMap,refIndex+16, Origin+Z*(height+wallThick),Z);  
  ModelSupport::buildPlane(SMap,refIndex+19, Origin+Z*(height+wallThick+voidThick),Z);  

  if (width>0) {
    ModelSupport::buildPlane(SMap, refIndex+111, Origin-X*(width/2.0), X);
    ModelSupport::buildPlane(SMap, refIndex+112, Origin+X*(width/2.0), X);
    ModelSupport::buildPlane(SMap, refIndex+121, Origin-X*(width/2.0+wallThick), X);
    ModelSupport::buildPlane(SMap, refIndex+122, Origin+X*(width/2.0+wallThick), X);
  }

  if (!(cellCooling)) return;

  double rRing; // ring radius
  int ringIndex(refIndex);
  for (int i=0; i<nRings; i++) {
    rRing = radius*(i+1)/(nRings+1);
    //    std::cout << "ring radius: " << rRing << std::endl;
    ModelSupport::buildCylinder(SMap, ringIndex+21, Origin, Z, rRing);
    ModelSupport::buildCylinder(SMap, ringIndex+22, Origin, Z, rRing+cellWallThick);
    ModelSupport::buildCylinder(SMap, ringIndex+23, Origin, Z, rRing-cellCoolThick);
    ModelSupport::buildCylinder(SMap, ringIndex+24, Origin, Z, rRing+cellWallThick+cellCoolThick);
    ringIndex += 10;
  }
  
  double hLayer = height/nLayers; // layer height
  double zLayer; // layer z-coordinate
  int layerIndex(refIndex);
  for (int i=0; i<nLayers-1; i++) {
    zLayer = hLayer*(i+1);
    std::cout << "layer: " << zLayer << std::endl;
    ModelSupport::buildPlane(SMap,layerIndex+25, Origin+Z*zLayer, Z);  
    ModelSupport::buildPlane(SMap,layerIndex+26, Origin+Z*(zLayer+cellWallThick), Z);  
    layerIndex += 10;
  }


  return; 
}

void
BeRef::addToInsertChain(attachSystem::ContainedComp& CC) const
  /*!
    Adds this object to the containedComp to be inserted.
    \param CC :: ContainedComp object to add to this
  */
{
  for(int i=refIndex+1;i<cellIndex;i++)
    CC.addInsertCell(i);
    
  return;
}


void
BeRef::createObjects(Simulation& System, const std::string &TargetSurfBoundary)
  /*!
    Create the vaned moderator
    \param System :: Simulation to add results
   */
{
  ELog::RegMethod RegA("BeRef","createObjects");

  std::string Out;
  
  // [2:1381] There are 2 types of cells: object cells (Monte Carlo objects = MC qhulls)

  if (cellCooling) {
    int ringIndex(refIndex);
    //    int layerIndex(refIndex);
    for (int i=0; i<=nRings; i++) {
      //      std::cout << "ring" << i << std::endl;
      if (i==(nRings)) {
	Out=ModelSupport::getComposite(SMap,refIndex, ringIndex, " -7 -6 15 14M");
	System.addCell(MonteCarlo::Qhull(cellIndex++,refMat,0.0,Out));
      } else if (i==0) {
	Out=ModelSupport::getComposite(SMap,refIndex, ringIndex, " -6 15 -23M");
	System.addCell(MonteCarlo::Qhull(cellIndex++,refMat,0.0,Out));
      } else {
	Out=ModelSupport::getComposite(SMap,refIndex, ringIndex, " -6 15 14M -23M");
	System.addCell(MonteCarlo::Qhull(cellIndex++,refMat,0.0,Out));
      }

      if (i!=(nRings)) {
	Out=ModelSupport::getComposite(SMap,refIndex, ringIndex, " -6 15 23M -21M");
	System.addCell(MonteCarlo::Qhull(cellIndex++, cellCoolMat, 0.0,Out));

	Out=ModelSupport::getComposite(SMap,refIndex, ringIndex, " -6 15 21M -22M");
	System.addCell(MonteCarlo::Qhull(cellIndex++, cellWallMat, 0.0,Out));

	Out=ModelSupport::getComposite(SMap,refIndex, ringIndex, " -6 15 22M -24M");
	System.addCell(MonteCarlo::Qhull(cellIndex++, cellCoolMat, 0.0,Out));
      }
      ringIndex += 10;
    }
  } else { // no cell cooling
    if (width<0) {
      Out=ModelSupport::getComposite(SMap,refIndex," -7 -6 15 ");
      System.addCell(MonteCarlo::Qhull(cellIndex++,refMat,0.0,Out));
    } else if (width>0) {
      Out=ModelSupport::getComposite(SMap,refIndex," -7 -6 15 111 -112 ");
      System.addCell(MonteCarlo::Qhull(cellIndex++,refMat,0.0,Out));
      Out=ModelSupport::getComposite(SMap,refIndex," -7 -6 15 -121 ");
      System.addCell(MonteCarlo::Qhull(cellIndex++, refMat1, 0.0, Out));
      Out=ModelSupport::getComposite(SMap,refIndex," -7 -6 15 122 ");
      System.addCell(MonteCarlo::Qhull(cellIndex++, refMat1, 0.0, Out));

      Out=ModelSupport::getComposite(SMap,refIndex," -7 -6 15 -111 +121 ");
      System.addCell(MonteCarlo::Qhull(cellIndex++, wallMat, 0.0, Out));
      Out=ModelSupport::getComposite(SMap,refIndex," -7 -6 15 -122 +112 ");
      System.addCell(MonteCarlo::Qhull(cellIndex++, wallMat, 0.0, Out));
    }
  }
  
  // reflector wall:
  Out=ModelSupport::getComposite(SMap,refIndex," -17 -16 (7:6:-15)"); // " -17 5 -16 (7:-5:6)"
  Out += TargetSurfBoundary;
  System.addCell(MonteCarlo::Qhull(cellIndex++,wallMat,0.0,Out));

  // clearance
  Out=ModelSupport::getComposite(SMap,refIndex," -18 -19 (17:16)");
  Out += TargetSurfBoundary;
  System.addCell(MonteCarlo::Qhull(cellIndex++, 0, 0.0, Out));
    
  // outer surface:
  Out=ModelSupport::getComposite(SMap,refIndex," -18 -19 ");
  Out += TargetSurfBoundary;
  addOuterSurf(Out); // mandatory that an object had an outer surface. it is also possible to add outer union. do not make it too big / complicated.
  
  return; 

}

void
BeRef::createLinks()
  /*!
    Creates a full attachment set
    Links/directions going outwards true.
  */
{
  /*
    Few lines of magick.
    [2:1680] How something else could connect to this object, or examine this object, or find out where it is, or someting.
    Defines how an object gets moved. Every object has to have them.
   */

  // [2:1850] The idea is that the link poin and the Z are the axis you want to connect along.
  // [2:1872] How to link/connect multiple surfacess (add...)

  FixedComp::setConnect(0, Origin+Y*(radius+wallThick+voidThick),-Y); // origin on the y-radius (i.e. somewhere on the beam) pointing inwards (because -Y). 99% of the links are pointing outwards. Actually, -Y or +Y does not matter so much, but the surface sign does matter. Stuart says it should be +Y here.
  FixedComp::setLinkSurf(0, SMap.realSurf(refIndex+18));

  FixedComp::setConnect(1, Origin+Y*(radius),-Y);
  FixedComp::setLinkSurf(1, SMap.realSurf(refIndex+7));

  //  FixedComp::setConnect(1,Origin-Z*(height/2.0+wallThick),-Z); // down in z and the link points in negative z. This is exactly how it should be.
  //  FixedComp::setLinkSurf(1,-SMap.realSurf(refIndex+15));

  FixedComp::setConnect(2, Origin+Z*(height/2.0+wallThick+voidThick),Z);
  FixedComp::setLinkSurf(2, SMap.realSurf(refIndex+19));

  FixedComp::setConnect(3, Origin+Y*(radius+wallThick),-Y);
  FixedComp::setLinkSurf(3, SMap.realSurf(refIndex+17));

  FixedComp::setConnect(4, Origin+Z*(height/2.0),Z);
  FixedComp::setLinkSurf(4, SMap.realSurf(refIndex+6));

  FixedComp::setConnect(5, Origin+Z*(wallThick),Z);
  FixedComp::setLinkSurf(5, SMap.realSurf(refIndex+15));


  return;
}


  void BeRef::createAll(Simulation& System, const attachSystem::FixedComp& FC, const attachSystem::FixedComp& TargetFC, const long int tIndex)
{
  /*!
    Extrenal build everything
    \param System :: Simulation
    \param FC :: FixedComponent for origin
    \param TargetFC :: target object
    \param tIndex :: link to upper/lower Target plane for upper/lower BeRef

    In our case the reflector has to be built relative to an origin and an axes set. 
    If you take a simple fixed object, then the axes is the axes set of this fixed object and the origin is the origin of this fixed object,
    which does not mean that the reflector and the object have the same origin, that's just the way you start and then you add the next bits.
   */

  ELog::RegMethod RegA("BeRef","createAll");
  // the order matters:

  populate(System.getDataBase()); // populate variables
  createUnitVector(FC); // take fixed component, then apply shift and angle rotation (transformation) for this object centre
  createSurfaces();

  const std::string TSurf=(tIndex<0) ? 
    TargetFC.getLinkComplement(static_cast<size_t>(-(tIndex+1))) : 
    TargetFC.getLinkString(static_cast<size_t>(tIndex));


  createObjects(System, TSurf);
  createLinks();
  insertObjects(System);       

  return;
}

}  // NAMESPACE instrumentSystem
