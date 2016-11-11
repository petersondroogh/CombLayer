/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   essBuildInc/Chicane.h
 *
 * Copyright (c) 2016 by Konstantin Batkov
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
#ifndef essSystem_Chicane_h
#define essSystem_Chicane_h

class Simulation;

namespace essSystem
{

/*!
  \class Chicane
  \version 1.0
  \author Konstantin Batkov
  \date DATE
  \brief Chicane
*/

class Chicane : public attachSystem::ContainedComp,
  public attachSystem::FixedOffset
{
 private:

  const std::string keyName;    ///< component name
  const int surfIndex;             ///< Index of surface offset
  int cellIndex;                ///< Cell index

  size_t nSegments;             ///< Number of segments
  std::vector<double> length;   ///< Lengths of all segments
  double width;                 ///< Width
  double height;                ///< Height

  int mat;                   ///< main material
  
  void populate(const FuncDataBase&);
  void createUnitVector(const attachSystem::FixedComp&);
  
  void createSurfaces(const attachSystem::FixedComp&,
		      const size_t&,const size_t&,const size_t&);
  void createLinks(const attachSystem::FixedComp&,
		   const size_t&,const size_t&,const size_t&);
  void createObjects(Simulation&,const attachSystem::FixedComp&,
		     const size_t&,const size_t&,const size_t&);

  
 public:

  Chicane(const std::string&,const int&);
  Chicane(const Chicane&);
  Chicane& operator=(const Chicane&);
  virtual ~Chicane();
  
  void createAll(Simulation&,const attachSystem::FixedComp&,
		 const attachSystem::FixedComp&,
		 const size_t&,const size_t&,const size_t&);

};

}

#endif
 

