 /********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   constructInc/SquareFMask.h
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
#ifndef xraySystem_SquareFMask_h
#define xraySystem_SquareFMask_h

class Simulation;

namespace xraySystem
{
  /*!
    \class SquareFMask
    \version 1.0
    \author S. Ansell
    \date June 2015
    \brief Variable detemine hole type
  */
  
class SquareFMask :
  public attachSystem::ContainedSpace,
  public attachSystem::FixedOffset,
  public attachSystem::CellMap,
  public attachSystem::SurfMap
{
 private:

  double width;                 ///< Main radius
  double height;                ///< Main height
  double length;                ///< thickness of collimator
  
  double innerAWidth;           ///< front width
  double innerAHeight;          ///< front height
  
  double minLength;               ///< point of min closure
  double innerMinWidth;           ///< min width at closure
  double innerMinHeight;          ///< min height at closure

  double innerBWidth;           ///< back width
  double innerBHeight;          ///< back height 

  double flangeAInRadius;        ///< Joining Flange inner radius
  double flangeAOutRadius;       ///< Joining Flange outer radius 
  double flangeALength;          ///< Joining Flange length

  double flangeBInRadius;        ///< Joining Flange inner radius
  double flangeBOutRadius;       ///< Joining Flange outer radius 
  double flangeBLength;          ///< Joining Flange length

  int flangeMat;                ///< material of flange
  int mat;                      ///< material
  int voidMat;                  ///< inner material
  
  void createUnitVector(const attachSystem::FixedComp&,
			const long int);
  void createSurfaces();
  void createObjects(Simulation&);
  void createLinks();
  
 public:
  
  SquareFMask(const std::string&);
  SquareFMask(const SquareFMask&);
  SquareFMask& operator=(const SquareFMask&);
  virtual ~SquareFMask() {}  ///< Destructor

  void populate(const FuncDataBase&);

  void createAll(Simulation&,const attachSystem::FixedComp&,
		 const long int);
  
};

}

#endif
 
