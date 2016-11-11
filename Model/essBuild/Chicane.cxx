/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   essBuild/Chicane.cxx
 *
 * Copyright (c) 2004-2016 by Konstantin Batkov
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
#include "FixedOffset.h"
#include "surfDBase.h"
#include "surfDIter.h"
#include "surfDivide.h"
#include "SurInter.h"
#include "mergeTemplate.h"

#include "Chicane.h"

namespace essSystem
{

  Chicane::Chicane(const std::string& Key, const int& N)  :
  attachSystem::ContainedComp(),
  attachSystem::FixedOffset(Key+StrFunc::makeString(N),6),
  keyName(Key+StrFunc::makeString(N)),
  surfIndex(ModelSupport::objectRegister::Instance().cell(keyName)),
  cellIndex(surfIndex+1)
  /*!
    Constructor BUT ALL variable are left unpopulated.
    \param Key :: Name for item in search
  */
{}

Chicane::Chicane(const Chicane& A) : 
  attachSystem::ContainedComp(A),
  attachSystem::FixedOffset(A),
  surfIndex(A.surfIndex),cellIndex(A.cellIndex),
  nSegments(A.nSegments),
  length(A.length),width(A.width),height(A.height),
  mat(A.mat)
  /*!
    Copy constructor
    \param A :: Chicane to copy
  */
{}

Chicane&
Chicane::operator=(const Chicane& A)
  /*!
    Assignment operator
    \param A :: Chicane to copy
    \return *this
  */
{
  if (this!=&A)
    {
      attachSystem::ContainedComp::operator=(A);
      attachSystem::FixedOffset::operator=(A);
      cellIndex=A.cellIndex;
      nSegments=A.nSegments;
      length=A.length;
      width=A.width;
      height=A.height;
      mat=A.mat;
    }
  return *this;
}

Chicane::~Chicane() 
  /*!
    Destructor
  */
{}

void
Chicane::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable data base
  */
{
  ELog::RegMethod RegA("Chicane","populate");

  FixedOffset::populate(Control);

  nSegments=Control.EvalVar<size_t>(keyName+"NSegments");
  for (size_t i=0; i<nSegments; i++)
    {
      const double l=Control.EvalVar<double>(StrFunc::makeString(keyName+"Length", i+1));
      length.push_back(l);
    }
  
  width=Control.EvalVar<double>(keyName+"Width");
  height=Control.EvalVar<double>(keyName+"Height");

  mat=ModelSupport::EvalMat<int>(Control,keyName+"Mat");

  return;
}
  
void
Chicane::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: object for origin
  */
{
  ELog::RegMethod RegA("Chicane","createUnitVector");

  FixedComp::createUnitVector(FC);
  applyShift(xStep,yStep,zStep);
  applyAngleRotate(xyAngle,zAngle);

  return;
}
  
void
Chicane::createSurfaces(const attachSystem::FixedComp& FC,
			const size_t& innerLP,
			const size_t& outerLP,
			const size_t& roofLP)
  /*!
    Create All the surfaces
    \param FC :: Bunker
    \param innerLP :: inner link ponit of Bunker
    \param outerLP :: outer link ponit of Bunker
    \param roofLP  :: link point to the inner roof of Bunker
  */
{
  ELog::RegMethod RegA("Chicane","createSurfaces");

  const Geometry::Plane *pRoof = SMap.realPtr<Geometry::Plane>(FC.getLinkSurf(roofLP));
  const Geometry::Cylinder *cylInner = SMap.realPtr<Geometry::Cylinder>(FC.getLinkSurf(innerLP));

  // bridge surface for createLinks:
  ModelSupport::buildPlane(SMap,surfIndex+1,Origin,Y);

  ModelSupport::buildCylinder(SMap,surfIndex+7,
			      cylInner->getCentre(),
			      cylInner->getNormal(),
			      cylInner->getRadius()-length[0]);

  int SI(surfIndex+10);
  double L(0.0);
  for (size_t i=1; i<=4; i++)
    {
      L += length[i];
      ModelSupport::buildCylinder(SMap,SI+7,
				  cylInner->getCentre(),
				  cylInner->getNormal(),
				  cylInner->getRadius()+L);
      SI += 10;
    }
  SMap.addMatch(surfIndex+57,FC.getLinkSurf(outerLP));

  ModelSupport::buildPlane(SMap,surfIndex+3,Origin-X*(width/2.0),X);
  ModelSupport::buildPlane(SMap,surfIndex+4,Origin+X*(width/2.0),X);

  SMap.addMatch(surfIndex+5,-FC.getLinkSurf(roofLP));
  ModelSupport::buildShiftedPlane(SMap,surfIndex+6,pRoof,height);
  ModelSupport::buildShiftedPlane(SMap,surfIndex+16,pRoof,height*2.0);
  
  return;
}
  
void
Chicane::createObjects(Simulation& System,
		       const attachSystem::FixedComp& FC,
		       const size_t& innerLP,
		       const size_t& outerLP,
		       const size_t& roofLP)
  /*!
    Adds the all the components
    \param System :: Simulation to create objects in
    \param FC :: Bunker
    \param innerLP :: inner link ponit of Bunker
    \param outerLP :: outer link ponit of Bunker
    \param roofLP  :: link point to the inner roof of Bunker
  */
{
  ELog::RegMethod RegA("Chicane","createObjects");

  const std::string innerSurf = FC.getLinkComplement(innerLP);
  const std::string outerSurf = FC.getLinkComplement(outerLP);
  const std::string roofSurf = FC.getLinkComplement(roofLP);

  std::string Out;
  const std::string side=ModelSupport::getComposite(SMap,surfIndex," 1 3 -4 ");
  
  Out=ModelSupport::getComposite(SMap,surfIndex," 5 -6 7 -27 ") + side;
  System.addCell(MonteCarlo::Qhull(cellIndex++,mat,0.0,Out));

  Out=ModelSupport::getComposite(SMap,surfIndex," 6 -16 17 -47 ") + side;
  System.addCell(MonteCarlo::Qhull(cellIndex++,mat,0.0,Out));

  Out=ModelSupport::getComposite(SMap,surfIndex," 5 -6 37 -57 ") + side;
  System.addCell(MonteCarlo::Qhull(cellIndex++,mat,0.0,Out));

  Out=ModelSupport::getComposite(SMap,surfIndex,
				 " ((5 -6 7 -27) : (6 -16 17 -47) : (5 -6 37 -57)) ") + side;
  addOuterSurf(Out);

  return;
}

  
void
Chicane::createLinks(const attachSystem::FixedComp& FC,
		     const size_t& innerLP,
		     const size_t& outerLP,
		     const size_t& roofLP)
  /*!
    Create all the links
    \param FC :: Bunker
    \param innerLP :: inner link ponit of Bunker
    \param outerLP :: outer link ponit of Bunker
    \param roofLP  :: link point to the inner roof of Bunker
  */
{
  ELog::RegMethod RegA("Chicane","createLinks");

  const Geometry::Plane *pRoof = SMap.realPtr
    <Geometry::Plane>(FC.getLinkSurf(roofLP));
  const Geometry::Cylinder *cylInner = SMap.realPtr
    <Geometry::Cylinder>(FC.getLinkSurf(innerLP));
  const Geometry::Cylinder *cylOuter = SMap.realPtr
    <Geometry::Cylinder>(FC.getLinkSurf(outerLP));

  FixedComp::setConnect(0,Origin+Z*pRoof->getDistance(),Z);
  FixedComp::setLinkSurf(0,-FC.getLinkSurf(roofLP));

  FixedComp::setConnect(1,Origin+Z*(pRoof->getDistance()+height*2.0),Z);
  FixedComp::setLinkSurf(1,SMap.realSurf(surfIndex+16));

  FixedComp::setConnect(2,Origin-X*(width/2.0),-X);
  FixedComp::setLinkSurf(2,-SMap.realSurf(surfIndex+3));

  FixedComp::setConnect(3,Origin+X*(width/2.0),X);
  FixedComp::setLinkSurf(3,SMap.realSurf(surfIndex+3));

  FixedComp::setConnect(4,Origin+Y*(cylInner->getRadius()-length[0]),-Y);
  FixedComp::setLinkSurf(4,-SMap.realSurf(surfIndex+7));
  FixedComp::setBridgeSurf(4,SMap.realSurf(surfIndex+1));

  FixedComp::setConnect(5,Origin+Y*(cylOuter->getRadius()),Y);
  FixedComp::setLinkSurf(5,SMap.realSurf(cylOuter->getName()));
  FixedComp::setBridgeSurf(5,SMap.realSurf(surfIndex+1));

  return;
}
  
  

  
void
Chicane::createAll(Simulation& System,
		   const attachSystem::FixedComp& origFC,
		   const attachSystem::FixedComp& FC,
		   const size_t& innerLP,
		   const size_t& outerLP,
		   const size_t& roofLP)
  /*!
    Generic function to create everything
    \param System :: Simulation item
    \param origFC :: Central origin
    \param FC :: Bunker
    \param innerLP :: inner link ponit of Bunker
    \param outerLP :: outer link ponit of Bunker
    \param roofLP  :: link point to the inner roof of Bunker
  */
{
  ELog::RegMethod RegA("Chicane","createAll");

  populate(System.getDataBase());
  createUnitVector(origFC);
  createSurfaces(FC, innerLP, outerLP, roofLP);
  createLinks(FC, innerLP, outerLP, roofLP);

  createObjects(System, FC, innerLP, outerLP, roofLP);
  insertObjects(System);              

  return;
}

}  // essSystem essSystem
