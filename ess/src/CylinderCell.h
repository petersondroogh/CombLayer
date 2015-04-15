#ifndef essSystem_CylinderCell_h
#define essSystem_CylinderCell_h

class Simulation;

namespace essSystem
{
  class CylinderCell : public attachSystem::ContainedComp, public attachSystem::FixedComp
    {
    private:

      const int refIndex;             ///< Index of surface offset
      int cellIndex;                  ///< Cell index. Every object has a space of about 10000 cells unless you request more.

      double xStep;                   ///< X step
      double yStep;                   ///< Y step
      double zStep;                   ///< Z step
      double xyAngle;                 ///< XY Angle
      double zAngle;                  ///< Z Angle

      double radius;                  ///< Radius
      double height;                  ///< Height
      double wallThick;               ///< Wall thickness
      double voidThick;               ///< Void clearance thickness

      int refMat;                     ///< reflector material - lower half
      //  int refMatTop;                     ///< reflector material - upper half
      int wallMat;                    ///< wall Material

      // Functions:

      void populate(const FuncDataBase&);
      void createUnitVector(const attachSystem::FixedComp&);

      void createSurfaces();
      void createObjects(Simulation&);
      void createLinks();

    public:

      CylinderCell(const std::string&);
      CylinderCell(const CylinderCell&);
      CylinderCell& operator=(const CylinderCell&);
      virtual ~CylinderCell();

      virtual void addToInsertChain(attachSystem::ContainedComp&) const;
      void createAll(Simulation&,const attachSystem::FixedComp&);

      inline int getRefIndex() { return refIndex; }
      inline int getCellIndex() { return cellIndex; }
  
    };

}

#endif
 
