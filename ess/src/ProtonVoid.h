#ifndef essSystem_ProtonVoid_h
#define essSystem_ProtonVoid_h

class Simulation;

namespace essSystem
{
  class ProtonVoid : public attachSystem::ContainedComp,
    public attachSystem::FixedComp
      {
      private:
  
	const int pvIndex;            ///< Index of surface offset
	int cellIndex;                ///< Cell index
	int protonVoidCell;           ///< Inner void cell

	double viewWidth;
	double viewHeight;
	double length1; // length of the first segment (near the target wheel)
	double radius2; // radius of the second segment (near the monolith side surface)
  
	void populate(const Simulation&);
	void createUnitVector(const attachSystem::FixedComp&);
	void createSurfaces();
	void createLinks();
	void createObjects(Simulation&,const std::string&,const std::string&);

      public:

	ProtonVoid(const std::string&);
	ProtonVoid(const ProtonVoid&);
	ProtonVoid& operator=(const ProtonVoid&);
	~ProtonVoid();

	int getVoidCell() const { return protonVoidCell; }
	void createAll(Simulation&,const attachSystem::FixedComp&,
		       const long int,
		       const attachSystem::FixedComp&,const long int);

	virtual void addToInsertChain(attachSystem::ContainedComp&) const; // has to be there [2:1040]
      };

}

#endif
