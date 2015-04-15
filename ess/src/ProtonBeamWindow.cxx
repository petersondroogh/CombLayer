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
#include "ProtonBeamWindow.h"


namespace essSystem
{

  ProtonBeamWindow::ProtonBeamWindow(const std::string& Key) :
    attachSystem::ContainedComp(),attachSystem::FixedComp(Key,6),
    pbwIndex(ModelSupport::objectRegister::Instance().cell(Key)),
    cellIndex(pbwIndex+1)
    /*!
      Constructor
      \param Key :: Name of construction key
    */
  {}


  ProtonBeamWindow::~ProtonBeamWindow()
  /*!
    Destructor
  */
  {}

  void ProtonBeamWindow::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable table to use
  */
  {
    ELog::RegMethod RegA("ProtonBeamWindow","populate");


    // Master values
    active=Control.EvalVar<int>(keyName+"Active");
    xStep=Control.EvalVar<double>(keyName+"XStep");
    yStep=Control.EvalVar<double>(keyName+"YStep");
    zStep=Control.EvalVar<double>(keyName+"ZStep");
    xyAngle=Control.EvalVar<double>(keyName+"XYangle");
    zAngle=Control.EvalVar<double>(keyName+"Zangle");

    linacMat=Control.EvalVar<int>(keyName+"LinacMat");
    monolithMat=Control.EvalVar<int>(keyName+"MonolithMat");


    pipeRadius=Control.EvalVar<double>(keyName+"PipeRadius"); // outer radius
    pipeWallThick=Control.EvalVar<double>(keyName+"PipeWallThick"); // along y
    pipeWallMat=Control.EvalVar<int>(keyName+"PipeWallMat");
    pipeCoolingMat=Control.EvalVar<int>(keyName+"PipeCoolingMat");

    frameHeight = Control.EvalVar<double>(keyName+"FrameHeight");
    frameWidth  = Control.EvalVar<double>(keyName+"FrameWidth");
    frameThick  = Control.EvalVar<double>(keyName+"FrameThick");
    frameMat    = Control.EvalVar<int>(keyName+"FrameMat");

    panpipeHeight = Control.EvalVar<double>(keyName+"PanpipeHeight");
    panpipeWidth  = Control.EvalVar<double>(keyName+"PanpipeWidth");

    nPipes = ceil(panpipeWidth/pipeRadius/2.0)+1; // add +1 in order to add a piece of an additional pipe in the case of shfit<0 in createObjects

    return;
  }

  void ProtonBeamWindow::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: Fixed Component
  */
  {
    ELog::RegMethod RegA("ProtonBeamWindow","createUnitVector");
    attachSystem::FixedComp::createUnitVector(FC);
    applyShift(xStep,yStep,zStep);
    applyAngleRotate(xyAngle,zAngle);
  
    return;
  }

  void ProtonBeamWindow::createSurfaces(const attachSystem::FixedComp& FC)
  {
    /*!
      Create Surfaces for the F5 collimator
    */
    ELog::RegMethod RegA("ProtonBeamWindow","createSurfaces");
      

    //    SMap.addMatch(pbwIndex+1, FC.getLinkSurf(2)); // inner radius
    // frame:
    ModelSupport::buildPlane(SMap,pbwIndex+1, Origin-X*(frameWidth/2.0), X);
    ModelSupport::buildPlane(SMap,pbwIndex+2, Origin+X*(frameWidth/2.0), X);
    ModelSupport::buildPlane(SMap,pbwIndex+3, Origin-Y*(frameThick/2.0), Y);
    ModelSupport::buildPlane(SMap,pbwIndex+4, Origin+Y*(frameThick/2.0), Y);
    ModelSupport::buildPlane(SMap,pbwIndex+5, Origin-Z*(frameHeight/2.0),Z);
    ModelSupport::buildPlane(SMap,pbwIndex+6, Origin+Z*(frameHeight/2.0),Z);

    ModelSupport::buildPlane(SMap,pbwIndex+10, Origin, Y); // a vertical plane through pipe centres to separate accelerator void and monolith He
    
    // panpipe section
    ModelSupport::buildPlane(SMap,pbwIndex+11, Origin-X*(panpipeWidth/2.0), X);
    ModelSupport::buildPlane(SMap,pbwIndex+12, Origin+X*(panpipeWidth/2.0), X);
    ModelSupport::buildPlane(SMap,pbwIndex+15, Origin-Z*(panpipeHeight/2.0), Z);
    ModelSupport::buildPlane(SMap,pbwIndex+16, Origin+Z*(panpipeHeight/2.0), Z);


    // pipes
    const double shift=-0.0001; // be sure the neighbouring pipes intersect a little bit in order to avoid geometry errors
    const double xmin = -panpipeWidth/2.0+pipeRadius+shift/2;
    for (int ipipe=0; ipipe<nPipes; ipipe++) {
      ModelSupport::buildCylinder(SMap, pbwIndex+100+ipipe, Origin+X*(ipipe*(pipeRadius*2.0+shift)+xmin), Z, pipeRadius);
      ModelSupport::buildCylinder(SMap, pbwIndex+200+ipipe, Origin+X*(ipipe*(pipeRadius*2.0+shift)+xmin), Z, pipeRadius-pipeWallThick);
    }

    return; 
  }

  void ProtonBeamWindow::createObjects(Simulation& System, const std::string &ProtonVoidBoundary)
  {
    ELog::RegMethod RegA("ProtonBeamWindow","createObjects");

    if (!active) return;

    std::string Out;
 
    //    int voidMat = 0;//1001;

    Out = ModelSupport::getComposite(SMap, pbwIndex, " 1 -2 3 -4 5 -6 ");
    addOuterSurf(Out);

    Out = ModelSupport::getComposite(SMap, pbwIndex, " 1 -2 3 -4 5 -6 (-11:12:-3:4:-15:16) "); // frame
    System.addCell(MonteCarlo::Qhull(cellIndex++, frameMat, 0.0, Out));

    // void from the accelerator side
    Out = ModelSupport::getComposite(SMap, pbwIndex, " 11 -12 3 -4 15 -16 -10");
    int pipeIndex(pbwIndex);
    for (int ipipe=0; ipipe<nPipes; ipipe++) {
      Out += ModelSupport::getComposite(SMap, pipeIndex, " +100 ");
      pipeIndex++;
    }
    //    Out += ProtonVoidBoundary;
    System.addCell(MonteCarlo::Qhull(cellIndex++, linacMat, 0.0, Out));

    // void from the monolith side
    Out = ModelSupport::getComposite(SMap, pbwIndex, " 11 -12 3 -4 15 -16 +10");
    pipeIndex = pbwIndex;
    for (int ipipe=0; ipipe<nPipes; ipipe++) {
      Out += ModelSupport::getComposite(SMap, pipeIndex, " +100 ");
      pipeIndex++;
    }
    //    Out += ProtonVoidBoundary;
    System.addCell(MonteCarlo::Qhull(cellIndex++, monolithMat, 0.0, Out));

    pipeIndex = pbwIndex;
    Out = ModelSupport::getComposite(SMap, pbwIndex, " 11 -12 15 -16 ");
    Out += "(";
    for (int ipipe=0; ipipe<nPipes; ipipe++) {
      Out += ModelSupport::getComposite(SMap, pipeIndex, " (-100 +200) ");
      if (ipipe!=nPipes-1) Out += " : ";
      pipeIndex++;
    }
    Out += ")";
    //Out += ProtonVoidBoundary;
    System.addCell(MonteCarlo::Qhull(cellIndex++, pipeWallMat, 0.0, Out));
      
    pipeIndex = pbwIndex;
    Out = ModelSupport::getComposite(SMap, pbwIndex, " 11 -12 15 -16 ");
    Out += "(";
    for (int ipipe=0; ipipe<nPipes; ipipe++) {
      Out += ModelSupport::getComposite(SMap, pipeIndex, " -200 ");
      if (ipipe!=nPipes-1) Out += " : ";
      pipeIndex++;
    }    
    Out += ")";
    //    Out += ProtonVoidBoundary;
    System.addCell(MonteCarlo::Qhull(cellIndex++, pipeCoolingMat, 0.0, Out));
    

    return; 
  }

  void ProtonBeamWindow::createLinks()
  /*!
    Creates a full attachment set
    Links/directions going outwards true.
  */
  {
    ELog::RegMethod RegA("ProtonBeamWindow","createLinks");

    FixedComp::setConnect(0,  Origin-Y*(frameThick/2), -Y);
    FixedComp::setLinkSurf(0, SMap.realSurf(pbwIndex+3));

    FixedComp::setConnect(0,  Origin+Y*(frameThick/2), Y);
    FixedComp::setLinkSurf(0, SMap.realSurf(pbwIndex+4));

    return;
  }


  void
  ProtonBeamWindow::createAll(Simulation& System,
			      const attachSystem::FixedComp& FC, const long int pvIndex)
  /*!
    pvIndex is the surface of proton void cylinder
  */
  {
    ELog::RegMethod RegA("ProtonBeamWindow","createAll");
    populate(System.getDataBase());

    createUnitVector(FC);
    createSurfaces(FC);

    const std::string pvSurf = (pvIndex<0) ? FC.getLinkComplement(static_cast<size_t>(-(pvIndex+1))) : FC.getLinkString(static_cast<size_t>(pvIndex));
    createObjects(System, pvSurf);
    createLinks();
    insertObjects(System);       

    return;
  }

  void ProtonBeamWindow::addToInsertChain(attachSystem::ContainedComp& CC) const
  /*!
    Adds this object to the containedComp to be inserted.
    \param CC :: ContainedComp object to add to this
  */
  {
    for(int i=pbwIndex+1;i<cellIndex;i++)
      CC.addInsertCell(i);
    
    return;
  }


}  // namespace essSystem
