/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   ESSBeam/odin/ODIN.h
 *
 * Copyright (c) 2004-2016 by Stuart Ansell
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
#ifndef essSystem_ODIN_h
#define essSystem_ODIN_h

namespace attachSystem
{
  class FixedComp;
  class FixedGroup;
  class TwinComp;
  class CellMap;
}

namespace constructSystem
{
  class Jaws;
  class DiskChopper;
  class ChopperPit;
  class RotaryCollimator;
  class PinHole;
  class ChopperUnit;
  class VacuumPipe;
}

namespace essSystem
{
  class Bunker;
  class BunkerInsert;
  class Hut;
  class RentrantBS;
  
  /*!
    \class ODIN
    \version 1.0
    \author S. Ansell
    \date January 2013
    \brief Main moderator system for ESS
  */
  
class ODIN : public attachSystem::CopiedComp
{
 private:

  /// Stop at [0:Complete / 1:Mono Wall / 2:Inner Bunker / 3:Outer Bunker ]
  int stopPoint;  

  /// Main Beam Axis [for construction]
  std::shared_ptr<attachSystem::FixedComp> odinAxis;

  /// Pipe between bunker and the wall
  std::shared_ptr<constructSystem::VacuumPipe> VPipeB;
  /// Elliptic guide from 5.5 to 6metre
  std::shared_ptr<beamlineSystem::GuideLine> FocusB;
  /// Quad chopper housing
  std::shared_ptr<constructSystem::ChopperUnit> ChopperA;
  /// Quad blades
  std::shared_ptr<constructSystem::DiskChopper> QDisk;

  /// Pipe after quad chopper
  std::shared_ptr<constructSystem::VacuumPipe> VPipeC;
  /// Basic tapered guid
  std::shared_ptr<beamlineSystem::GuideLine> FocusC;

  
  /// Tapper Unit
  std::shared_ptr<beamlineSystem::GuideLine> GuideA;
  /// T0 chopper [9-9.5m]
  std::shared_ptr<constructSystem::DiskChopper> T0Chopper;
  /// Tapper Unit
  std::shared_ptr<beamlineSystem::GuideLine> GuideB;
  /// Bunker insert
  std::shared_ptr<essSystem::BunkerInsert> BInsert;
  /// Guide in the Bunker wall
  std::shared_ptr<beamlineSystem::GuideLine> GuideC;
  /// Guide after the Bunker to first chopper
  std::shared_ptr<beamlineSystem::GuideLine> GuideD;

  /// Chopper pit for first outer bunker chopper
  std::shared_ptr<constructSystem::ChopperPit> PitA;
  /// Guide to Chopper to exterior
  std::shared_ptr<beamlineSystem::GuideLine> GuidePitAFront;
  /// Guide from Chopper to exterior
  std::shared_ptr<beamlineSystem::GuideLine> GuidePitABack;

  /// Guide from Chopper A to exterior
  std::shared_ptr<beamlineSystem::GuideLine> GuideE;


  /// Chopper pit for second choppers:
  std::shared_ptr<constructSystem::ChopperPit> PitB;
  /// Guide from Chopper to exterior [target]
  std::shared_ptr<beamlineSystem::GuideLine> GuidePitBFront;
  /// Guide from Chopper to exterior [Hutch side]
  std::shared_ptr<beamlineSystem::GuideLine> GuidePitBBack;
  /// Guide from Chopper to exterior
  std::shared_ptr<constructSystem::DiskChopper> ChopperB;
  /// Guide from chopper B to exterior
  std::shared_ptr<beamlineSystem::GuideLine> GuideF;


  /// Chopper pit for third choppers:
  std::shared_ptr<constructSystem::ChopperPit> PitC;
  /// Guide from Chopper to exterior [target]
  std::shared_ptr<beamlineSystem::GuideLine> GuidePitCFront;
  /// Guide from Chopper to exterior [Hutch side]
  std::shared_ptr<beamlineSystem::GuideLine> GuidePitCBack;


  /// Guide from chopper C to exterior
  std::shared_ptr<beamlineSystem::GuideLine> GuideG;

  /// The Hutch
  std::shared_ptr<Hut> Cave;
  /// Guide in Hutch 
  std::shared_ptr<beamlineSystem::GuideLine> GuideH;
  
  /// Collimator
  std::shared_ptr<constructSystem::PinHole> PinA;

  /// BeamStop
  std::shared_ptr<RentrantBS> BeamStop;

  void setBeamAxis(const attachSystem::FixedGroup&,const bool);
  
 public:
  
  ODIN(const std::string&);
  ODIN(const ODIN&);
  ODIN& operator=(const ODIN&);
  ~ODIN();
  
  void build(Simulation&,const attachSystem::FixedGroup&,
	     const Bunker&,const int);

};

}

#endif
