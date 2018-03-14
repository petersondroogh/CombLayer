/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   geometry/surfImplicates.cxx
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
#include <vector>
#include <map>
#include <list>
#include <stack>
#include <string>
#include <algorithm>

#include "Exception.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Triple.h"
#include "Quaternion.h"
#include "support.h"
#include "regexBuild.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "Surface.h"
#include "Quadratic.h"
#include "ArbPoly.h"
#include "CylCan.h"
#include "Cylinder.h"
#include "Cone.h"
#include "EllipticCyl.h"
#include "General.h"
#include "MBrect.h"
#include "Plane.h"
#include "Sphere.h"
#include "Torus.h"
#include "surfEqual.h"
#include "surfImplicates.h"


namespace Geometry
{

surfImplicates::surfImplicates() 
  //  functionMap({"PlanePlane",&surfImplicates::planePlane})
  /*!
    Constructor
  */
{
  functionMap.emplace("PlanePlane",&surfImplicates::planePlane);
}

surfImplicates& 
surfImplicates::Instance() 
  /*!
    Effective this object			
    \return surfImplicates object
  */
{
  static surfImplicates A;
  return A;
}

int
surfImplicates::isImplicate(const Surface* ASPtr,
			    const Surface* BSPtr) const
  /*!
    Determine if the pair is identical
    \param ASPtr :: first surface
    \param BSPtr :: second surface
    \return -1 / 1 if implicate 0 if not
  */
{
  ELog::RegMethod RegA("surfaceImplicates","isImplicate");

  const std::string AType=ASPtr->className();
  const std::string BType=BSPtr->className();

  const std::string combineA=AType+BType;
  STYPE::const_iterator mc=functionMap.find(combineA);
  if (mc!=functionMap.end())
    {
      testImp IPtr=mc->second;
      return (this->*IPtr)(ASPtr,BSPtr);
    }

  // all else failed
  return 0;
}

int
surfImplicates::planePlane(const Geometry::Surface* APtr,
			   const Geometry::Surface* BPtr) const
  /*!
    Determine if two planes are implicates
    \param APtr :: First plane pointer
    \param BPtr :: second plane pointer
   */
{
  // we already know these are planes
  const Plane* APlane=dynamic_cast<const Plane*>(APtr);
  const Plane* BPlane=dynamic_cast<const Plane*>(BPtr);

  const Geometry::Vec3D& ANorm=APlane->getNormal();
  const Geometry::Vec3D& BNorm=BPlane->getNormal();
  const double AD=APlane->getDistance();
  const double BD=BPlane->getDistance();
  
  if (ANorm==BNorm)
    return AD>BD ? -1 : 1;
  
  if (ANorm==-BNorm)
    return AD > BD ? 1 : -1;

  return 0;

}
  
} // NAMESPACE Geometry
