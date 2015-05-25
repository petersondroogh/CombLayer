#ifndef essSystem_BeRef_h
#define essSystem_BeRef_h

#include "CellMap.h"

class Simulation;

namespace essSystem
{

  /*!
    \class BeRef
    \author S. Ansell
    \version 1.0
    \date February 2013
    \brief Reflector object 
  */

  class BeRef : public attachSystem::ContainedComp, public attachSystem::FixedComp, public attachSystem::CellMap
    {
    private:

      const int refIndex;             ///< Index of surface offset
      int cellIndex;                  ///< Cell index. Every object has a space of about 10000 cells unless you request more.

      double xStep;                   ///< X step
      double yStep;                   ///< Y step
      double zStep;                   ///< Z step
      double xyAngle;                 ///< XY Angle
      double zAngle;                  ///< Z Angle

      double radius;                  ///< total Radius
      double height;                  ///< Height
      double wallThick;               ///< Wall thickness
      double voidThick;               ///< Void clearance thickness

      int refMat;                     ///< reflector material
      int wallMat;                    ///< wall Material

      double innerRadius;               ///< radius of the inner reflector part
      double innerHeight;               ///< height of the inner reflector part
      double innerWallThick;          ///< inner wall thickness
      int    innerRefMat;             ///< material of the inner reflector part (if i  
      int    innerWallMat;              ///< inner wall material


      double width;                   ///< reflector width in the case we wish to limit it by planes. Has no effect if < 0
      double length;                  ///< reflector length in the case we wish to limit it by planes. Has no effect if < 0
      int    refMat1;                 ///< material inside the reflector cylinder and outside width/length - e.g. in the area complementing refMat. used only if width>0 || length>0

      double VoidCellHeight;             ///< height (thickness) of void layer between target and reflector (WaterDisc)
      int VoidCellMat;             ///< material of (void by default) layer between target and reflector (WaterDisc)

      // Functions:

      void populate(const FuncDataBase&);
      void createUnitVector(const attachSystem::FixedComp&);

      void createSurfaces();
      void createObjects(Simulation&, const std::string&);
      void createLinks();

    public:

      BeRef(const std::string&);
      BeRef(const BeRef&);
      BeRef& operator=(const BeRef&);
      virtual ~BeRef();

      //      int getMainCell() const { return refIndex+1; } // Be
      //      int getWallCell() const { return refIndex+2; } // wall
      virtual void addToInsertChain(attachSystem::ContainedComp&) const; // has to be there [2:1040]
      void createAll(Simulation&,const attachSystem::FixedComp&, const attachSystem::FixedComp&, const long int); // [2:1070]

      inline int getRefIndex() { return refIndex; }
      inline int getCellIndex() { return cellIndex; }
  
    };

}

#endif
 
