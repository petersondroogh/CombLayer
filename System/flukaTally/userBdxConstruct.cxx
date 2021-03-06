/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   flukaTally/userBdxConstruct.cxx
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
#include "surfRegister.h"
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
#include "BaseMap.h"
#include "SurfMap.h"
#include "CellMap.h"
#include "LinkSupport.h"
#include "inputParam.h"

#include "Object.h"
#include "SimFLUKA.h"
#include "particleConv.h"
#include "flukaGenParticle.h"
#include "TallySelector.h"
#include "flukaTally.h"
#include "userBdx.h"
#include "userBdxConstruct.h" 


namespace flukaSystem
{

bool
userBdxConstruct::checkLinkCells(const Simulation& System,
				 const int cellA,
				 const int cellB)
  /*!
    Determine if two cells are connected via a common surface
    \param System :: Simulation for objects
    \param cellA :: Cell A
    \param cellB :: Cell B
  */
{
  ELog::RegMethod RegA("userBdxConstruct","checkLinkCells");
  const MonteCarlo::Object* APtr=System.findObject(cellA);
  const MonteCarlo::Object* BPtr=System.findObject(cellB);
  if (!APtr || !BPtr)
    return 0;

  const std::set<int>& ASurf=APtr->getSurfSet();
  const std::set<int>& BSurf=APtr->getSurfSet();

  for(const int CN : ASurf)
    if (BSurf.find(-CN)!=BSurf.end())
      return 1;

  return 0;
}

bool
userBdxConstruct::constructLinkRegion(const Simulation& System,
				      const std::string& FCname,
				      const std::string& FCindex,
				      int& cellA,int& cellB)
  /*!
    Construct a link region exiting the FixedComp link unit
    \param System :: Simulation to use	
    \param FCname :: name of fixed comp
    \param FCiindex :: name of link point
  */
{
  ELog::RegMethod RegA("userBdxConstruct","constructLinkRegion");
  
  const attachSystem::FixedComp* FCPtr=
    System.getObject<attachSystem::FixedComp>(FCname);

  if (!FCPtr) return 0;

  if (!FCPtr->hasSideIndex(FCindex)) return 0;
  const long int FCI=FCPtr->getSideIndex(FCindex);

  const int surfN=FCPtr->getLinkSurf(FCI);
  if (!surfN) return 0;

  const std::pair<const MonteCarlo::Object*,
	    const MonteCarlo::Object*> RefPair=
    System.findCellPair(FCPtr->getLinkPt(FCI),surfN);
  
  if (RefPair.first && RefPair.second)
    {
      cellA=RefPair.first->getName();
      cellB=RefPair.second->getName();
      return 1;
    }
  return 0;
}

bool
userBdxConstruct::constructSurfRegion(const Simulation& System,
				      const std::string& FCname,
				      const std::string& surfName,
				      const size_t indexA,
				      const size_t indexB,
				      int& cellA,int& cellB)
  /*!
    Construct a link region exiting the SurfMap link unit
    FCname also names a groupRange which is used 
    to ensure that cellA is part of the groupRange
    \param System :: Simulation to use	
    \param FCname :: name of SurfMap
    \param surfName :: name of surface [signed]
    \param indexA :: Index of region found in primary
    \param indexB :: Index region found in secondary
    \param cellA :: Primary region cell number
    \param cellB :: Secondary region cell number
  */
{
  ELog::RegMethod RegA("userBdxConstruct","constructSurfRegion");
  
  const attachSystem::SurfMap* SMPtr=
    System.getObject<attachSystem::SurfMap>(FCname);

  if (!SMPtr || surfName.empty()) return 0;
  
  const int surfN=SMPtr->getSignedSurf(surfName);
  if (!surfN) return 0;
  // throws on error [unlikely because SurfMap is good]
  const groupRange& activeGrp=System.getGroup(FCname);
  const std::pair<const MonteCarlo::Object*,
	    const MonteCarlo::Object*> RefPair=
    System.findCellPair(surfN,activeGrp,indexA,indexB);
  
  if (RefPair.first && RefPair.second)
    {
      cellA=RefPair.first->getName();
      cellB=RefPair.second->getName();
      return 1;
    }
  return 0;
}

void 
userBdxConstruct::createTally(SimFLUKA& System,
			      const std::string& PType,const int fortranTape,
			      const int cellA,const int cellB,
			      const bool eLog,const double Emin,
			      const double Emax,const size_t nE,
			      const bool aLog,const double Amin,
			      const double Amax,const size_t nA)
  /*!
    An amalgamation of values to determine what sort of mesh to put
    in the system.
    \param System :: SimFLUKA to add tallies
    \param fortranTape :: output stream
    \param CellA :: initial region
    \param CellB :: secondary region
    \param eLog :: energy in log bins
    \param aLog :: angle in log bins
    \param Emin :: Min energy 
    \param Emax :: Max energy 
  */
{
  ELog::RegMethod RegA("userBdxConstruct","createTally");

  const flukaGenParticle& FG=flukaGenParticle::Instance();
    
  userBdx UD(fortranTape);
  UD.setParticle(FG.nameToFLUKA(PType));

  UD.setCell(cellA,cellB);
  UD.setEnergy(eLog,Emin,Emax,nE);
  UD.setAngle(aLog,Amin,Amax,nA);
  
  System.addTally(UD);

  return;
}


void
userBdxConstruct::processBDX(SimFLUKA& System,
			     const mainSystem::inputParam& IParam,
			     const size_t Index) 
  /*!
    Add BDX tally (s) as needed
    - Input:
    -- particle FixedComp index
    -- particle cellA  cellB
    -- particle SurfMap name
    \param System :: SimFLUKA to add tallies
    \param IParam :: Main input parameters
    \param Index :: index of the -T card
  */
{
  ELog::RegMethod RegA("userBdxConstruct","processBdx");

  
  const std::string particleType=
    IParam.getValueError<std::string>("tally",Index,1,"tally:ParticleType");
  const std::string FCname=
    IParam.getValueError<std::string>("tally",Index,2,"tally:Object/Cell");
  const std::string FCindex=
    IParam.getValueError<std::string>("tally",Index,3,"tally:linkPt/Cell");

  size_t itemIndex(4);
  int cellA(0);
  int cellB(0);
  if (
      (!StrFunc::convert(FCname,cellA) ||
       !StrFunc::convert(FCindex,cellB) ||
       !checkLinkCells(System,cellA,cellB) ) &&
      
      !constructLinkRegion(System,FCname,FCindex,cellA,cellB)
      )
    
    {
      // special class because must give regions
      itemIndex+=2;
      const size_t regionIndexA=IParam.getDefValue(0,"tally",Index,4);
      const size_t regionIndexB=IParam.getDefValue(0,"tally",Index,5);

      if (!constructSurfRegion(System,FCname,FCindex,
			       regionIndexA,regionIndexB,cellA,cellB))
	throw ColErr::InContainerError<std::string>
	  (FCname+":"+FCindex,"No connecting surface on regions");
    }
  
  ELog::EM<<"Regions connected from "<<cellA<<" "<<cellB<<ELog::endDiag;  

  // This needs to be more sophisticated
  const int nextId=System.getNextFTape();
  
  const double EA=IParam.getDefValue<double>(1e-9,"tally",Index,itemIndex++);
  const double EB=IParam.getDefValue<double>(1000,"tally",Index,itemIndex++);
  const size_t NE=IParam.getDefValue<size_t>(200,"tally",Index,itemIndex++); 

  const double AA=IParam.getDefValue<double>(0.0,"tally",Index,itemIndex++);
  const double AB=IParam.getDefValue<double>(2*M_PI,"tally",Index,itemIndex++);
  const size_t NA=IParam.getDefValue<size_t>(1,"tally",Index,itemIndex++); 

  
  userBdxConstruct::createTally(System,particleType,nextId,
				cellA,cellB,
				1,EA,EB,NE,
				0,AA,AB,NA);
  
  return;      
}  
  
void
userBdxConstruct::writeHelp(std::ostream& OX) 
  /*!
    Write out help
    \param OX :: Output stream
  */
{
  OX<<
    "recordType filename \n"
    "  --recordType avaiable:\n"
    "    source : trajectory : local : continuous\n"
    "    sourceLoss : trajLoss : user";

  return;
}

}  // NAMESPACE flukaSystem
