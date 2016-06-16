/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   photonInc/RecColl.h
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
#ifndef photonSystem_RecColl_h
#define photonSystem_RecColl_h

class Simulation;

namespace photonSystem
{

struct LInfo
{
  size_t nDisk;                ///< number of units
  double thick;                ///< Thickness
  std::vector<double> Radii;   ///< Radii
  std::vector<int> Mat;        ///< Material
  std::vector<double> Temp;    ///< Temperature
  
  void resize(const size_t);

};

/*!
  \class RecColl
  \author S. Ansell
  \version 1.0
  \date Janurary 2015
  \brief Specialized for a layered cylinder pre-mod object
*/

class RecColl : public attachSystem::ContainedComp,
  public attachSystem::FixedOffset,
  public attachSystem::CellMap
{
 private:

  const int layerIndex;         ///< Index of surface offset
  int cellIndex;                ///< Cell index


  double outerRadius;                ///< Outer radius
  size_t nLayers;                    ///< Layer count
  std::vector<LInfo> LVec;           ///< Layer Info

  void populate(const FuncDataBase&);
  void createUnitVector(const attachSystem::FixedComp&,
			const long int);

  void createSurfaces();
  void createObjects(Simulation&);
  void createLinks();
  
 public:

  RecColl(const std::string&);
  RecColl(const RecColl&);
  RecColl& operator=(const RecColl&);
  virtual ~RecColl();
  virtual RecColl* clone() const;
  
  void createAll(Simulation&,const attachSystem::FixedComp&,
		 const long int);
};

}

#endif
 
