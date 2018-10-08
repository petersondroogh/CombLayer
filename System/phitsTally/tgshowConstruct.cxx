/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   phitsTally/tgshowConstruct.cxx
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
#include <iterator>
#include <memory>

#include "Exception.h"
#include "FileReport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "GTKreport.h"
#include "OutputLog.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "support.h"
#include "stringCombine.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "Rules.h"
#include "HeadRule.h"
#include "Code.h"
#include "varList.h"
#include "FuncDataBase.h"
#include "groupRange.h"
#include "objectGroups.h"
#include "Simulation.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "LinkSupport.h"
#include "inputParam.h"

#include "SimPHITS.h"
#include "particleConv.h"
#include "TallySelector.h"
#include "meshConstruct.h"
#include "MeshXYZ.h"
#include "phitsTally.h"
#include "TGShow.h"
#include "tgshowConstruct.h" 



namespace phitsSystem
{

void 
tgshowConstruct::createTally(SimPHITS& System,
			     const int ID,
			     const Geometry::Vec3D& APt,
			     const Geometry::Vec3D& BPt,
			     const std::array<size_t,3>& MPts) 
/*!
    An amalgamation of values to determine what sort of mesh to put
    in the system.
    \param System :: SimPHITS to add tallies
    \param ID :: output stream
    \param APt :: Lower point 
    \param BPt :: Upper point 
    \param MPts :: Points ot use
  */
{
  ELog::RegMethod RegA("tgshowConstruct","createTally");

  TGShow UB(ID);
  //  UB.setParticle(PType);
  UB.setCoordinates(APt,BPt);
  UB.setIndex(MPts);
  //  System.addTally(UB);

  return;
}

std::string
tgshowConstruct::convertTallyType(const std::string& TType)
  /*!
    Construct to convert words into an index
    \param TType :: incoming value / particle
    \return particle name [upper case] or string of number
  */
{
  ELog::RegMethod RegA("tgshowConstruct","convertTallyType");
  const particleConv& pConv=particleConv::Instance();
  
  if (TType=="help")
    {
      ELog::EM<<"Tally Type:"<<ELog::endDiag;
      ELog::EM<<" -- Particle"<<ELog::endDiag;
      throw ColErr::ExitAbort("Help termianted");
    }
      
  
  std::ostringstream cx;
  
  if (pConv.hasName(TType))
    return StrFunc::toUpperString(pConv.nameToPHITS(TType));

  throw ColErr::InContainerError<std::string>(TType,"TType not in TMap");
}

void
tgshowConstruct::processMesh(SimPHITS& System,
			      const mainSystem::inputParam& IParam,
			      const size_t Index) 
  /*!
    Add mesh tally (s) as needed
    \param System :: SimPHITS to add tallies
    \param IParam :: Main input parameters
    \param Index :: index of the -T card
  */
{
  ELog::RegMethod RegA("tgshowConstruct","processMesh");
  
  const std::string tallyType=
    IParam.getValueError<std::string>("tally",Index,1,"tallyType");
  const std::string tallyParticle=
    tgshowConstruct::convertTallyType(tallyType);

  // This needs to be more sophisticated
  const int nextId=10; // System.getNextFTape();
  
  const std::string PType=
    IParam.getValueError<std::string>("tally",Index,2,"object/free"); 
  Geometry::Vec3D APt,BPt;
  std::array<size_t,3> Nxyz;
  
  if (PType=="object")
    tallySystem::meshConstruct::getObjectMesh
      (System,IParam,"tally",Index,3,APt,BPt,Nxyz);
  else if (PType=="free")
    tallySystem::meshConstruct::getFreeMesh
      (IParam,"tally",Index,3,APt,BPt,Nxyz);

  tgshowConstruct::createTally(System,nextId,APt,BPt,Nxyz);
  
  return;      
}
  
void
tgshowConstruct::writeHelp(std::ostream& OX) 
  /*!
    Write out help
    \param OX :: Output stream
  */
{
  OX<<
    "tallyType free Vec3D Vec3D Nx Ny Nz \n"
    "tallyType object objectName LinkPt  Vec3D Vec3D Nx Ny Nz \n"
    "  -- Object-link point is used to construct basis set \n"
    "     Then the Vec3D are used as the offset points ";

  return;
}

}  // NAMESPACE tallySystem
