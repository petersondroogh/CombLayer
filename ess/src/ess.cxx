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
#include <array>
#include <algorithm>
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>

#include "Exception.h"
#include "MersenneTwister.h"
#include "FileReport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "GTKreport.h"
#include "OutputLog.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "InputControl.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Tensor.h"
#include "Vec3D.h"
#include "inputParam.h"
#include "Triple.h"
#include "NRange.h"
#include "NList.h"
#include "Tally.h"
#include "TallyCreate.h"
#include "Transform.h"
#include "Quaternion.h"
#include "localRotate.h"
#include "masterRotate.h"
#include "Surface.h"
#include "Quadratic.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Line.h"
#include "Rules.h"
#include "surfIndex.h"
#include "Code.h"
#include "varList.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "ModeCard.h"
#include "PhysCard.h"
#include "LSwitchCard.h"
#include "PhysImp.h"
#include "KGroup.h"
#include "SrcData.h"
#include "Source.h"
#include "KCode.h"
#include "PhysicsCards.h"
#include "BasicWWE.h"
#include "MainProcess.h"
#include "SimProcess.h"
#include "SimInput.h"
#include "SurInter.h"
#include "Simulation.h"
#include "SimPHITS.h"
#include "PointWeights.h"
#include "ContainedComp.h"
#include "ContainedGroup.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "LinearComp.h"
#include "mainJobs.h"
#include "Volumes.h"
#include "DefPhysics.h"
#include "variableSetup.h"
#include "ImportControl.h"
#include "SourceCreate.h"
#include "SourceSelector.h"
#include "TallySelector.h"
#include "World.h"
#include "makeESS.h"

MTRand RNG(12345UL);

///\cond STATIC
namespace ELog 
{
  ELog::OutputLog<EReport> EM;
  ELog::OutputLog<FileReport> FM("Spectrum.log");
  ELog::OutputLog<FileReport> RN("Renumber.txt");   ///< Renumber
  ELog::OutputLog<StreamReport> CellM;
}
///\endcond STATIC

int 
main(int argc,char* argv[])
{
  int exitFlag(0);                // Value on exit
  ELog::RegMethod RControl("","main");
  mainSystem::activateLogging(RControl);

  std::string Oname;
  std::vector<std::string> Names;  
  std::map<std::string,std::string> Values;  
  std::map<std::string,std::string> AddValues;  
  std::map<std::string,double> IterVal;           // Variable to iterate 

  // PROCESS INPUT:
  InputControl::mainVector(argc,argv,Names);
  mainSystem::inputParam IParam;
  createESSInputs(IParam);

  // IParam.regDefItem<std::string>("lowMod","lowMod",1,std::string("lowMod"));
  // IParam.setDesc("lowMod","Type of low modeartor");

  // IParam.regDefItem<std::string>("topMod","topMod",1,std::string("topMod"));
  // IParam.setDesc("topMod","Type of top modeartor");

  IParam.regDefItem<std::string>("topModFlowGuide","topModFlowGuide",1,std::string("topModFlowGuide"));
  IParam.setDesc("topModFlowGuide","Type of top moderator cooling (Onion|Standard|None)");

  IParam.regDefItem<std::string>("lowModFlowGuide","lowModFlowGuide",1,std::string("lowModFlowGuide"));
  IParam.setDesc("lowModFlowGuide","Type of low moderator cooling (Onion|Standard|None)");

  IParam.regDefItem<std::string>("topPreCooling","topPreCooling",1,std::string("topPreCooling"));
  IParam.setDesc("topPreCooling","Type of top premoderator cooling (Onion|None)");

  IParam.regDefItem<std::string>("lowWaterDisc","lowWaterDisc",1,std::string("lowWaterDisc"));
  IParam.setDesc("lowWaterDisc","Water disc betwen the wheel and the bottom moderator (On|Off)");

  IParam.regDefItem<std::string>("topWaterDisc","topWaterDisc",1,std::string("topWaterDisc"));
  IParam.setDesc("topWaterDisc","Water disc betwen the wheel and the top moderator (On|Off)");

  IParam.regDefItem<std::string>("thermalCylMod","thermalCylMod",1,std::string("thermalCylMod"));
  IParam.setDesc("thermalCylMod","Thermal moderator of cylinder shape betwen the wheel/waterCylinder and the bottom [cold] moderator (On|Off)");

  IParam.regDefItem<int>("nF5", "nF5", 1,2);
  IParam.setDesc("nF5","Number of F5 collimators to build. The collimators will be named as F5, F15, F25 etc. The corresponding variables must exist.");

  IParam.regFlag("rotate", "rotate");
  IParam.setDesc("rotate","Rotate to the Alan's coordinate system");

  IParam.regFlag("matmesh", "matmesh");
  IParam.setDesc("matmesh","Generate material mesh (to compare two geometries)");
  
  const int iteractive(IterVal.empty() ? 0 : 1);   
  Simulation* SimPtr=createSimulation(IParam,Names,Oname);
  if (!SimPtr) return -1;

  // The big variable setting
  setVariable::EssVariables(SimPtr->getDataBase());
  InputModifications(SimPtr,IParam,Names);
  
  // Definitions section 
  int MCIndex(0);
  const int multi=IParam.getValue<int>("multi");
  try
    {
      while(MCIndex<multi)
	{
	  if (MCIndex)
	    {
	      ELog::EM.setActive(4);    // write error only
	      ELog::FM.setActive(4);    
	      ELog::RN.setActive(0);    
	      // if (iteractive)
	      // 	mainSystem::incRunTimeVariable
	      // 	  (SimPtr->getDataBase(),IterVal);
	    }

	  SimPtr->resetAll();

	  essSystem::makeESS ESSObj;
	  World::createOuterObjects(*SimPtr);
	  ESSObj.build(SimPtr,IParam);

	  SDef::sourceSelection(*SimPtr,IParam);
	  // kbat start
	  physicsSystem::PhysicsCards &PC=SimPtr->getPC();
		  /*SDef::Source A = PC.getSDefCard();
	  //	  A.setComp("y", -590);
	  
	  SDef::SrcData D3(3);
	  SDef::SrcData D4(4);
	  
	  SDef::SrcInfo SI3('H');
	  SI3.addData(-7.0);
	  SI3.addData(7.0);
	  
	  SDef::SrcProb SP3(1);
	  SP3.addData(0.0);
	  SP3.addData(1.0);

	  D3.addUnit(SI3);
	  D3.addUnit(SP3);


	  SDef::SrcInfo SI4('H');
	  SI4.addData(-1.6);
	  SI4.addData(1.6);

	  SDef::SrcProb SP4(4);
	  SP4.addData(0.0);
	  SP4.addData(1.0);

	  D4.addUnit(SI4);
	  D4.addUnit(SP4);
	  
	  A.setData("x", D3);
	  A.setData("z", D4);
	  A.setActive();*/
	  //	  A.write(std::cout);



	  SimPtr->removeComplements();
	  SimPtr->removeDeadSurfaces(0);         

	  ModelSupport::setDefaultPhysics(*SimPtr,IParam);
	  // kbat start
	  //	  PC.setMode("n p | h / z");  // why does not work?
	  PC.setPrdmp("2j 01 02");

	  physicsSystem::LSwitchCard &lea=PC.getLEA();
	  lea.setValues("lca", "2 1 1 23 1 1 0 1 0 1"); // 9J 01
	  //	  lea.setValues("lcb", "3490.0 3490.0 2490.0 2490.0 800.0 800.0 -1.0 -1.0"); // 8J
	  //	  lea.setValues("lcc", "1.0 45.0"); // 2J
	  lea.setValues("lea", "1 4 1 0 1 0 0 1"); // 8J

	  //	  SimPtr->getPC().setPrintNum("010  040  070  098  102  110  130  140  170"); // does not work
	  physicsSystem::PhysCard &nCard = PC.addPhysCard("cut", "n");
	  double neutronCutoff = 0.0;
	  nCard.setValues(4, 1e+8, neutronCutoff, -0.5, -0.25); // same as Alan's J 0.0, but Stuart sets 4,1e+8, 0, 0.4, -0.1

	  physicsSystem::PhysCard &pCard = PC.addPhysCard("cut", "p");
	  pCard.setValues(2, neutronCutoff, 0.0);

	  physicsSystem::PhysCard &allCard = PC.addPhysCard("cut", "h / d t s a | z");
	  allCard.setValues(2, neutronCutoff, 0.0);

	  // kbat end
	  const int renumCellWork=tallySelection(*SimPtr,IParam);
	  // kbat start
	  if (IParam.flag("rotate")) { // rotate to the Alan's coordinate system
	    // SA says I need to use master rotation here - see hist t2 model
	    masterRotate& MR = masterRotate::Instance();
	    MR.addRotation(Geometry::Vec3D(1,0,0), Geometry::Vec3D(0,0,0), -90.0);
	    MR.addRotation(Geometry::Vec3D(0,1,0), Geometry::Vec3D(0,0,0), 180.0);
	    //	  MR.addMirror(Geometry::Plane(1,0,Geometry::Vec3D(0,0,0), Geometry::Vec3D(1,0,0)));
	  }
	  // kbat end
	  SimPtr->masterRotation();
	  if (createVTK(IParam,SimPtr,Oname))
	    {
	      delete SimPtr;
	      ModelSupport::objectRegister::Instance().reset();
	      return 0;
	    }
	  if (IParam.flag("endf"))
	    SimPtr->setENDF7();
	  //	  createMeshTally(IParam,SimPtr);

	  // nps - number of incident particles (guessed by kbat)
	  SimPtr->getPC().setNPS(20000);
	  //	  SimPtr->getPC().setPrintNum("10 11");

	  SimProcess::importanceSim(*SimPtr,IParam);
	  SimProcess::inputPatternSim(*SimPtr,IParam); // energy cut etc

	  if (renumCellWork)
	    tallyRenumberWork(*SimPtr,IParam);
	  tallyModification(*SimPtr,IParam);

	  if (IParam.flag("cinder"))
	    SimPtr->setForCinder();

	  // // Cut energy tallies:
	  // if (IParam.flag("ECut"))
	  //   SimPtr->setEnergy(IParam.getValue<double>("ECut"));

	  // Ensure we done loop
	  do
	    {
	      SimProcess::writeIndexSim(*SimPtr,Oname,MCIndex);
	      MCIndex++;
	    }
	  while(!iteractive && MCIndex<multi);
	}
      exitFlag=SimProcess::processExitChecks(*SimPtr,IParam);
      if (IParam.flag("cinder"))
	SimPtr->writeCinder();
      ModelSupport::calcVolumes(SimPtr,IParam);
      ModelSupport::objectRegister::Instance().write("ObjectRegister.txt");
    }
  catch (ColErr::ExitAbort& EA)
    {
      if (!EA.pathFlag())
	ELog::EM<<"Exiting from "<<EA.what()<<ELog::endCrit;
      exitFlag=-2;
    }
  catch (ColErr::ExBase& A)
    {
      ELog::EM<<"EXCEPTION FAILURE :: "
	      <<A.what()<<ELog::endCrit;
      exitFlag= -1;
    }
  delete SimPtr;
  ModelSupport::objectRegister::Instance().reset();
  ModelSupport::surfIndex::Instance().reset();
  return exitFlag;
}
