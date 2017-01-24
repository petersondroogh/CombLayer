/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   commonVarInc/TwinGenerator.h
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
#ifndef essSystem_TwinGenerator_h
#define essSystem_TwinGenerator_h

class Simulation;

namespace setVariable
{

/*!
  \class TwinGenerator
  \version 1.0
  \author S. Ansell
  \date May 2016
  \brief TwinGenerator for variables
*/

class TwinGenerator
{
 private:

  double height;           ///< Full height
  double width;            ///< Full width
  double mainRadius;       ///< full chopper void radius
  double windowThick;      ///< window material thickness
  double ringRadius;       ///< Radius of the rinng 

  double motorRadius;      ///< Motor radius
  double motorOuter;       ///< Motor outer
  double portRadius;       ///< Port radius
  double portOuter;        ///< Port outer

  double portWidth;        ///< Port width
  double portHeight;       ///< Port height
  double portBoltStep;     ///< Port outer
  
  std::string wallMat;     ///< Main body material
  std::string portMat;     ///< Port material
  std::string sealMat;     ///< seal material
  std::string windowMat;    ///< seal material
  
 public:

  TwinGenerator();
  TwinGenerator(const TwinGenerator&);
  TwinGenerator& operator=(const TwinGenerator&);
  ~TwinGenerator();


  void setMaterial(const std::string&,const std::string&);
  void setFrame(const double,const double);  
  void setMainRadius(const double);


  void generateChopper(FuncDataBase&,
		       const std::string&,
		       const double,const double,
		       const double);
};

}

#endif
 