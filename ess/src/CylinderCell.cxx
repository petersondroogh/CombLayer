/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   essBuild/CylinderCell.cxx
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
#include "CylinderCell.h"


namespace essSystem
{

CylinderCell::CylinderCell(const std::string& Key) :
  attachSystem::ContainedComp(),attachSystem::FixedComp(Key,6),
  refIndex(ModelSupport::objectRegister::Instance().cell(Key)),
  cellIndex(refIndex+1)
  /*!
    Constructor
    \param Key :: Name of construction key
  */
{}


CylinderCell::~CylinderCell()
  /*!
    Destructor
   */
{}

void
CylinderCell::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable table to use
  */
{
  ELog::RegMethod RegA("CylinderCell","populate");


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

  return;
}

void
CylinderCell::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: Fixed Component
  */
{
  ELog::RegMethod RegA("CylinderCell","createUnitVector");
  attachSystem::FixedComp::createUnitVector(FC);
  applyShift(xStep,yStep,zStep);
  applyAngleRotate(xyAngle,zAngle);
  
  return;
}

void
CylinderCell::createSurfaces()
  /*!
    Create Surfaces for the Be
  */
{
  ELog::RegMethod RegA("CylinderCell","createSurfaces");
      
  ModelSupport::buildCylinder(SMap,refIndex+1, Origin,Z,radius);   // SMap - local surface map [2:1300]
  ModelSupport::buildPlane(SMap,refIndex+2, Origin-Z*(height/2), Z);  
  ModelSupport::buildPlane(SMap,refIndex+3, Origin+Z*(height/2), Z); // upper/lower plane of bottom Al container (other plane is the target upper/lower surface)
  //  ModelSupport::buildCylinder(SMap,refIndex+11, Origin,Z,radius+wallThick);
  //  ModelSupport::buildPlane(SMap,refIndex+12, Origin-Z*(height/2-wallThick), Z);
  //  ModelSupport::buildPlane(SMap,refIndex+13 Origin+Z*(height/2+wallThick), Z);

  return; 
}

void
CylinderCell::addToInsertChain(attachSystem::ContainedComp& CC) const
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
CylinderCell::createObjects(Simulation& System)
  /*!
    Create the vaned moderator
    \param System :: Simulation to add results
   */
{
  ELog::RegMethod RegA("CylinderCell","createObjects");

  std::string Out;
  
  // [2:1381] There are 2 types of cells: object cells (Monte Carlo objects = MC qhulls)

  Out=ModelSupport::getComposite(SMap,refIndex," -1 -3 2");
  System.addCell(MonteCarlo::Qhull(cellIndex++,refMat,0.0,Out));
    
  addOuterSurf(Out); // mandatory that an object had an outer surface. it is also possible to add outer union. do not make it too big / complicated.
  
  return; 

}

void
CylinderCell::createLinks()
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


  void CylinderCell::createAll(Simulation& System, const attachSystem::FixedComp& FC)
{
  /*!
    Extrenal build everything
    \param System :: Simulation
    \param FC :: FixedComponent for origin
   */

  ELog::RegMethod RegA("CylinderCell","createAll");
  // the order matters:

  populate(System.getDataBase()); // populate variables
  createUnitVector(FC); // take fixed component, then apply shift and angle rotation (transformation) for this object centre
  createSurfaces();

  createObjects(System);
  createLinks();
  insertObjects(System);       

  return;
}

}  // NAMESPACE instrumentSystem
