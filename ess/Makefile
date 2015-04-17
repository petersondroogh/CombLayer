#CL=../

INC=-I$(CL)/supportInc -I$(CL)/tallyInc -I$(CL)/logInc -I$(CL)/include -I$(CL)/inputInc -I$(CL)/globalInc -I$(CL)/monteInc -I$(CL)/processInc -I$(CL)/worldInc -I$(CL)/attachCompInc -I$(CL)/funcBaseInc -I$(CL)/geomInc -I$(CL)/buildInc -I$(CL)/chipInc -I$(CL)/physicsInc -I$(CL)/constructInc -I$(CL)/mersenneInc -I$(CL)/weightsInc -I$(CL)/sourceInc -I./src

CXXFLAGS= -std=c++11 -g -fPIC -Wconversion -W -Wall -Wextra -Wno-comment -fexceptions

# obj/ess.o should not be in OBJ:
OBJ=obj/essVariables.o obj/BulkModule.o obj/BeRef.o obj/CylinderCell.o obj/essDBMaterial.o obj/ConicModerator.o  obj/WedgeItem.o obj/BlockAddition.o obj/CylModerator.o  obj/GuideBay.o obj/SegWheel.o obj/WheelBase.o  obj/CylPreMod.o obj/essMod.o  obj/GuideItem.o obj/ShutterBay.o obj/Wheel.o obj/F5Collimator.o obj/FlightLine.o obj/ProtonVoid.o obj/RecModerator.o obj/ButterflyModerator.o obj/ProtonBeamWindow.o obj/OnionCooling.o obj/F5Calc.o obj/WaterPipe.o obj/Grooving.o obj/ThermalModerator.o obj/makeESS.o

CC=g++

all: ess

obj/essVariables.o: src/essVariables.cxx
	@echo "Compiling $@"
	@$(CC) -c $(CXXFLAGS) $(INC) $< -o $@

obj/ess.o: src/ess.cxx
	@echo "Compiling $@"
	$(CC) -c $(CXXFLAGS) $(INC) $< -o $@

obj/%.o: src/%.cxx src/%.h
	@echo "Compiling $@"
	@$(CC) -c $(CXXFLAGS) $(INC) $< -o $@

libEssBuild.so: $(OBJ)
	@echo "Linking $@"
	@$(CC) -std=c++11 -fPIC -shared -g -o libEssBuild.so $(OBJ)

ess: libEssBuild.so obj/ess.o
	@echo "Linking $@"
	@$(CC) $(CXXFLAGS) -o $@ $^ -fPIC -Wl,-rpath,$(CL) $(CL)/lib/libVisit.so   $(CL)/lib/libSrc.so $(CL)/lib/libPhysics.so   $(CL)/lib/libInput.so $(CL)/lib/libSource.so   $(CL)/lib/libMonte.so $(CL)/lib/libFuncBase.so   $(CL)/lib/libLog.so $(CL)/lib/libTally.so   $(CL)/lib/libConstruct.so   $(CL)/lib/libCrystal.so   $(CL)/lib/libTransport.so   $(CL)/lib/libScatMat.so $(CL)/lib/libMd5.so   $(CL)/lib/libEndf.so $(CL)/lib/libProcess.so   $(CL)/lib/libWorld.so $(CL)/lib/libWork.so   $(CL)/lib/libMonte.so $(CL)/lib/libGeometry.so   $(CL)/lib/libMersenne.so $(CL)/lib/libSrc.so   $(CL)/lib/libXml.so $(CL)/lib/libPoly.so   $(CL)/lib/libSupport.so $(CL)/lib/libWeights.so   $(CL)/lib/libGlobal.so   $(CL)/lib/libAttachComp.so   $(CL)/lib/libVisit.so $(CL)/lib/libSimMC.so -lgsl -lgslcblas -lm  -lboost_regex 


clean:
	@rm -fv obj/*.o libEssBuild.so ess

cleandata:
	@rm -fv  ObjectRegister.txt Spectrum.log Renumber.txt foo.* *.x runtp* out? *.ps *.pdf comou? *.ps *.root *~


butterfly: ess
	pstudy-nowrap -i xml/butterfly-static.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -fv flat1.x
	./ess -r -lowMod Butterfly -topMod Butterfly -nF5 6  -topModCooling None  -lowWaterDisc On -topWaterDisc On  -topPreCooling None -x flat.xml -T surface object Bulk 3 -TMod energy 1 "5E-9 20E-9"  butterfly
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' butterfly1.x
	sed -i -e "s/fq1   f d u s m e t c/*c1:n 90 5 0/" butterfly1.x
#	remove-weights butterfly1.x

butterfly-pstudy: ess
	pstudy-nowrap -i xml/butterfly.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -fv flat1.x
	./ess -r -lowMod Butterfly -topMod Butterfly -nF5 6  -topModCooling None  -lowWaterDisc On -topWaterDisc On  -topPreCooling None -x flat.xml -T surface object Bulk 3 -TMod energy 1 "5E-9 20E-9"  butterfly
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' butterfly1.x
	sed -i -e "s/fq1   f d u s m e t c/*c1:n 90 5 0/" butterfly1.x
#	remove-weights butterfly1.x




test: ess testButterfly
	@rm -f tdr1.x
	@./ess -r -x xml/tdr.xml --validCheck 10000  tdr > /dev/null && echo "TDR test passed"
	@rm -f pipe3011.x
	@pstudy-nowrap -i xml/flat-master-PIPE301.xml -setup > /dev/null
	@/bin/mv -f case001/inp pipe301.xml && rm -fr case001 > /dev/null
	@./ess -r -lowMod None -topModCooling Onion -x pipe301.xml  --validCheck 10000 pipe301 > /dev/null && echo "PIPE301 test passed"
	@rm -f flat1.x
	@pstudy-nowrap -i xml/flat-master.xml -setup > /dev/null
	@/bin/mv -f case001/inp flat.xml && rm -fr case001
	@./ess -r -topModCooling Onion -x flat.xml  --validCheck 10000 flat > /dev/null && echo "flat test passed"
	@pstudy-nowrap -i xml/flat-master.xml -setup > /dev/null
	@/bin/mv -f case001/inp flatpipeTube.xml  && rm -fr case001 > /dev/null
	@rm -f flatpipeTube1.x
	@./ess -r -lowMod Tube -v LowBeRefWallMat 13000 -v LowBeRefRefMat 4000 -v LowBeRefVoidThick 1 -topModCooling Onion -topPreCooling Onion -x flatpipeTube.xml --validCheck 10000 flatpipeTube > /dev/null && echo "flatpipeTube test passed"
	@rm -f flatpipeTube1.x

testButterfly: ess
	@pstudy-nowrap -i xml/butterfly.xml -setup > /dev/null
	@/bin/mv -f case001/inp flat.xml  && rm -fr case001 > /dev/null
	@./ess -r -nF5 6 -topMod M3A -topModCooling Onion -lowMod Butterfly -waterCylinder On  -topPreCooling None -x flat.xml -v LowButterflyFlightLineType 0 -v LowButterflyTopPreType 0 --validCheck 10000 butterfly > /dev/null && echo "Butterfly00 test passed"
	@./ess -r -nF5 6 -topMod M3A -topModCooling Onion -lowMod Butterfly -waterCylinder On  -topPreCooling None -x flat.xml -v LowButterflyFlightLineType 0 -v LowButterflyTopPreType 1 --validCheck 10000 butterfly > /dev/null && echo "Butterfly01 test passed"
	@./ess -r -nF5 6 -topMod M3A -topModCooling Onion -lowMod Butterfly -waterCylinder On  -topPreCooling None -x flat.xml -v LowButterflyFlightLineType 0 -v LowButterflyTopPreType 2 --validCheck 10000  butterfly > /dev/null && echo "Butterfly02 test passed"
	@./ess -r -nF5 6 -topMod M3A -topModCooling Onion -lowMod Butterfly -waterCylinder On  -topPreCooling None -x flat.xml -v LowButterflyFlightLineType 1 -v LowButterflyTopPreType 0 --validCheck 10000  butterfly > /dev/null && echo "Butterfly10 test passed"
	@./ess -r -nF5 6 -topMod M3A -topModCooling Onion -lowMod Butterfly -waterCylinder On  -topPreCooling None -x flat.xml -v LowButterflyFlightLineType 1 -v LowButterflyTopPreType 1 --validCheck 10000  butterfly > /dev/null && echo "Butterfly11 test passed"
	@./ess -r -nF5 6 -topMod M3A -topModCooling Onion -lowMod Butterfly -waterCylinder On  -topPreCooling None -x flat.xml -v LowButterflyFlightLineType 1 -v LowButterflyTopPreType 2 --validCheck 10000  butterfly > /dev/null && echo "Butterfly12 test passed"
	@rm -f butterfly1.x



tags:
	@etags src/*.cxx $(CL)/*/*.cxx

tar:
	@tar jcf ess-batkov.tar.bz2 src/*.h src/*.cxx Makefile xml ~/bin/pstudy-nowrap
