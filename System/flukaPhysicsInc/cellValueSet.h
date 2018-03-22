/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   physicsInc/cellValueSet.h
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
#ifndef flukaSystem_cellValueSet_h
#define flukaSystem_cellValueSet_h

namespace flukaSystem
{
  
/*!
  \class cellValueSet
  \version 1.0
  \date March 2018
  \author S.Ansell
  \brief Processes the physics cards in the FLUKA output
*/

template <size_t N>
class cellValueSet 
{
 private:

  /// Data type
  typedef  std::array<double,N> valTYPE;
  /// map type
  typedef  std::map<int,valTYPE> dataTYPE;
 
  const std::string keyName;               ///< Key name
  const std::string outName;               ///< Output name for FLUKA
  double whatValue;                        ///< What [1] value 
  
  dataTYPE dataMap;   ///< Values for cell

  bool cellSplit(const std::vector<int>&,
		 std::vector<std::tuple<int,int>>&,
		 std::vector<valTYPE>&) const;
  
 public:

  cellValueSet(const std::string&,const std::string&);
  cellValueSet(const std::string&,const std::string&,const double);
  cellValueSet(const cellValueSet&);
  cellValueSet& operator=(const cellValueSet&);
  virtual ~cellValueSet();

  void clearAll();

  void setValue(const int,const size_t,const double);
  void setValues(const int,const double);    
  void setValues(const int,const double,const double);
  void setValues(const int,const double,const double,const double);
  void writeFLUKA(std::ostream&,const std::vector<int>&,
		  const std::string&) const;

};

}

#endif
