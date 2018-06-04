/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   essBuild/DTLArray.cxx
 *
 * Copyright (c) 2018 by Konstantin Batkov
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
#include "surfDBase.h"
#include "surfDIter.h"
#include "surfDivide.h"
#include "SurInter.h"
#include "mergeTemplate.h"
#include "CellMap.h"
#include "CopiedComp.h"

#include "DTL.h"
#include "DTLArray.h"

namespace essSystem
{

DTLArray::DTLArray(const std::string& baseKey,const std::string& Key)  :
  attachSystem::CopiedComp(Key,Key),
  attachSystem::FixedOffset(newName,6),
  baseName(baseKey)
  /*!
    Constructor BUT ALL variable are left unpopulated.
    \param Key :: Name for item in search
  */
{}

DTLArray::DTLArray(const DTLArray& A) : 
  attachSystem::CopiedComp(A),
  attachSystem::FixedOffset(A),
  baseName(A.baseName),
  nDTL(A.nDTL),
  dtl(A.dtl)
  /*!
    Copy constructor
    \param A :: DTLArray to copy
  */
{}

DTLArray&
DTLArray::operator=(const DTLArray& A)
  /*!
    Assignment operator
    \param A :: DTLArray to copy
    \return *this
  */
{
  if (this!=&A)
    {
      nDTL=A.nDTL;
      dtl=A.dtl;
    }
  return *this;
}

DTLArray*
DTLArray::clone() const
/*!
  Clone self
  \return new (this)
 */
{
    return new DTLArray(*this);
}
  
DTLArray::~DTLArray() 
  /*!
    Destructor
  */
{}

void
DTLArray::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable data base
  */
{
  ELog::RegMethod RegA("DTLArray","populate");

  FixedOffset::populate(Control);
  nDTL=Control.EvalPair<int>(baseName,"","NDTLTanks");
  ELog::EM << nDTL << ELog::endDiag;

  return;
}
  
void
DTLArray::createAll(Simulation& System,
		       const attachSystem::FixedComp& FC,
		       const long int sideIndex)
  /*!
    Generic function to create everything
    \param System :: Simulation item
    \param FC :: Central origin
    \param sideIndex :: link point for origin
  */
{
  ELog::RegMethod RegA("DTLArray","createAll");

  populate(System.getDataBase());

  return;
}

}  // essSystem
