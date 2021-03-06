/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   bibBuildInc/WaterMod.h
 *
 * Copyright (c) 2004-2017 by Stuart Ansell
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
#ifndef bibSystem_bibWaterMod_h
#define bibSystem_bibWaterMod_h

class Simulation;

namespace bibSystem
{

/*!
  \class WaterMod
  \author S. Ansell
  \version 1.0
  \date April 2013
  \brief Specialized for wheel
*/

class WaterMod : public attachSystem::ContainedComp,
    public attachSystem::FixedOffset
{
 private:
  

  double width;                   ///< width of moderator
  double height;                 
 ///< height of moderator
  double depth;                   ///< depth of moderator
  double wallThick;               ///< Wall thickness

  double sideGap;                 ///< Gaps
  double frontGap;                ///< Gaps
  double backGap;                 ///< Gaps
  double vertGap;                 ///< Gaps

  double waterTemp;               ///< water temperature

  int waterMat;                   ///< Water material
  int wallMat;                    ///< Wall material

  // Functions:

  void populate(const FuncDataBase&);
  void createUnitVector(const attachSystem::FixedComp&,const long int);

  void createSurfaces();
  void createObjects(Simulation&,const attachSystem::ContainedComp&);
  void createLinks();

  public:

  WaterMod(const std::string&);
  WaterMod(const WaterMod&);
  WaterMod& operator=(const WaterMod&);
  virtual ~WaterMod();

  void createAll(Simulation&,const attachSystem::FixedComp&,
		 const long int,const attachSystem::ContainedComp&);
  
};

}

#endif
 
