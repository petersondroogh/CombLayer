#CL=$(HOME)/misc/soft/comblayer/ansell3/CombLayer
#CL=$(HOME)/misc/soft/comblayer/CombLayer
CL=$(CLSYS)

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

install: libEssBuild.so
	ln -s $(PWD)/libEssBuild.so ~/usr/local/lib/

clean:
	@rm -fv obj/*.o libEssBuild.so ess

cleandata:
	@rm -fv  ObjectRegister.txt Spectrum.log Renumber.txt foo.* *.x runtp* out? *.ps *.pdf comou? *.ps *.root *~

changelog:
	@svn2cl

view: ess
	rm -fv /tmp/essview*x
	./ess -r --validCheck 1000 /tmp/essview
	perl -pi -e "s/mt/c mt/g" /tmp/essview1.x
	rm -f foo.*; 
	test -e zoom && mcnpx name=foo. i=/tmp/essview1.x ip com=zoom
	test -e zoom || mcnpx name=foo. i=/tmp/essview1.x ip

github: ess
	/bin/cp -f src/*.h ~/misc/soft/comblayer/kbat/essBuildInc
	/bin/cp -f src/*.cxx ~/misc/soft/comblayer/kbat/essBuild
	/bin/mv -f ~/misc/soft/comblayer/kbat/essBuild/ess.cxx ~/misc/soft/comblayer/kbat/Main/ess.cxx

inp: ess
	@./ess -r -T point free  -353.553 353.553 18 -TMod energy 5 "2.04511e-10 2.08651e-10 2.12919e-10 2.17319e-10 2.21856e-10 2.26538e-10 2.31369e-10 2.36356e-10 2.41506e-10 2.46827e-10 2.52325e-10 2.58009e-10 2.63887e-10 2.69969e-10 2.76263e-10 2.8278e-10 2.8953e-10 2.96525e-10 3.03776e-10 3.11297e-10 3.191e-10 3.27201e-10 3.35614e-10 3.44356e-10 3.53444e-10 3.62896e-10 3.72733e-10 3.82975e-10 3.93645e-10 4.04767e-10 4.16368e-10 4.28474e-10 4.41116e-10 4.54326e-10 4.68139e-10 4.82591e-10 4.97722e-10 5.13577e-10 5.30202e-10 5.47647e-10 5.65967e-10 5.85223e-10 6.05478e-10 6.26803e-10 6.49275e-10 6.72977e-10 6.98002e-10 7.24448e-10 7.52427e-10 7.82058e-10 8.13475e-10 8.46824e-10 8.82266e-10 9.1998e-10 9.60166e-10 1.00304e-09 1.04886e-09 1.09789e-09 1.15043e-09 1.20684e-09 1.26751e-09 1.33287e-09 1.40341e-09 1.47971e-09 1.5624e-09 1.65223e-09 1.75002e-09 1.85677e-09 1.97359e-09 2.10179e-09 2.2429e-09 2.39871e-09 2.57135e-09 2.76331e-09 2.9776e-09 3.21782e-09 3.48832e-09 3.79442e-09 4.14266e-09 4.54114e-09 5e-09 5.53209e-09 6.15386e-09 6.88668e-09 7.75863e-09 8.8073e-09 1.00839e-08 1.16596e-08 1.36357e-08 1.61603e-08 2e-08 2.38755e-08 2.99899e-08 3.87932e-08 5.21288e-08 7.37394e-08 1e-07 1.91032e-07 3.95231e-07 1.25303e-06 2.60204e-05" -w 2.0 -TW a
	@perl -pi -e "s/^t5/c t5/" a1.x
#	perl -pi -e "s/nps 10000/stop f5 1E-2/" a1.x

onemod: ess
	@./ess -r -lowMod None -v BeRefRefMatLow 26316 -T point free  -353.553 353.553 18 -TMod energy 5 "1e-12 2.04511e-10 2.08651e-10 2.12919e-10 2.17319e-10 2.21856e-10 2.26538e-10 2.31369e-10 2.36356e-10 2.41506e-10 2.46827e-10 2.52325e-10 2.58009e-10 2.63887e-10 2.69969e-10 2.76263e-10 2.8278e-10 2.8953e-10 2.96525e-10 3.03776e-10 3.11297e-10 3.191e-10 3.27201e-10 3.35614e-10 3.44356e-10 3.53444e-10 3.62896e-10 3.72733e-10 3.82975e-10 3.93645e-10 4.04767e-10 4.16368e-10 4.28474e-10 4.41116e-10 4.54326e-10 4.68139e-10 4.82591e-10 4.97722e-10 5.13577e-10 5.30202e-10 5.47647e-10 5.65967e-10 5.85223e-10 6.05478e-10 6.26803e-10 6.49275e-10 6.72977e-10 6.98002e-10 7.24448e-10 7.52427e-10 7.82058e-10 8.13475e-10 8.46824e-10 8.82266e-10 9.1998e-10 9.60166e-10 1.00304e-09 1.04886e-09 1.09789e-09 1.15043e-09 1.20684e-09 1.26751e-09 1.33287e-09 1.40341e-09 1.47971e-09 1.5624e-09 1.65223e-09 1.75002e-09 1.85677e-09 1.97359e-09 2.10179e-09 2.2429e-09 2.39871e-09 2.57135e-09 2.76331e-09 2.9776e-09 3.21782e-09 3.48832e-09 3.79442e-09 4.14266e-09 4.54114e-09 5e-09 5.53209e-09 6.15386e-09 6.88668e-09 7.75863e-09 8.8073e-09 1.00839e-08 1.16596e-08 1.36357e-08 1.61603e-08 2e-08 2.38755e-08 2.99899e-08 3.87932e-08 5.21288e-08 7.37394e-08 1e-07 1.91032e-07 3.95231e-07 1.25303e-06 2.60204e-05" -TMod time 5 "1 2 3" -w 2.0 -TW onemod

d2mod: ess
	@./ess -r  -lowMod Tube -v BeRefRefMatLow 82000 -v BulkMat2 82000 -T point free  -196 -1.26 -25.0 -v F5ColXStep -196 -v F5ColYStep -1.26 -v F5ColZStep -25.0 -v F5ColLength 100 -v F5ColXB -12.39 -v F5ColYB -12.33 -v F5ColZB -14.7 -v F5ColXC -12.46 -v F5ColYC 12.50 -v F5ColZC -14.7 -v F5ColZG -35.3 -TMod energy 5 "1.854971E-08 2.266040E-08" -TMod time 5 "0 200i 2E+6" -w 2.0 -TW d2mod


tdr: ess
	rm -fv tdr1.x
	./ess -r -x xml/tdr.xml --validCheck 1000  tdr
	rm -f foo.*; 
	test -e zoom && mcnpx name=foo. i=tdr1.x ip com=zoom
	test -e zoom || mcnpx name=foo. i=tdr1.x ip

onlyTarget: ess
	@./ess -r -v BeRefRefMatLow 0 -v BeRefRefMatTop 0 -v BeRefWallMat 0 -v BulkMat1 0 -v BulkMat2 0 -v BulkMat3 0 -v LowAFlightLinerMat1 0 -v LowBFlightLinerMat1 0 -v LowModMat 0 -v LowModMaterial1 0 -v LowModMaterial2 0 -v LowModMaterial3 0 -v LowModMaterial4 0 -v LowModMaterial5 0 -v LowModMaterial6 0 -v LowPreABlockWallMat 0 -v LowPreABlockWaterMat 0 -v LowPreBBlockWallMat 0 -v LowPreBBlockWaterMat 0 -v LowPreMaterial1 0 -v LowPreMaterial2 0 -v LowPreMaterial3 0 -v LowPreMaterial4 0 -v ShutterBayMat 0 -v TopAFlightLinerMat1 0 -v TopModMat 0 -v TopModMaterial1 0 -v TopModMaterial2 0 -v TopModMaterial3 0 -v TopModMaterial4 0 -v TopModMaterial5 0 -v TopModMaterial6 0 -v TopPreABlockWallMat 0 -v TopPreABlockWaterMat 0 -v TopPreBBlockWallMat 0 -v TopPreBBlockWaterMat 0 -v TopPreMaterial1 0 -v TopPreMaterial2 0 -v TopPreMaterial3 0 -v TopPreMaterial4 0 -v WheelCladShaftMat 0 -v WheelMainShaftMat 0 -v ProtonBeamWindowActive 0  -T surface object TopMod 1 -TMod energy 0 "1e-11 251log 3000" onlyTarget
	@perl -pi -e "s/e1    1e-11 251log 3000/`cat tallies.inp`/" onlyTarget1.x


pipe301: ess
	pstudy-nowrap -i xml/flat-master-PIPE301.xml -setup
	/bin/mv -f case001/inp pipe301.xml && rm -fr case001
	@rm -f pipe3011.x
	./ess -r -lowMod None -topModCooling Onion -x pipe301.xml  -v WheelMainShaftMat 26000 -v WheelCladShaftMat 26000  -T point free  489.016592713104 0 14.45 -TMod energy 5 "5E-9 2E-8 1E-7 4E-7  1.0 3E+3" -TMod time 5 "1E+8"  --validCheck 1000 pipe301
# -v F5ColLength 288.9 -v F5ColXStep 489.016592713104 -v F5ColYStep 0 -v F5ColZStep 14.45 -v F5ColXB 11.6397 -v F5ColYB -5.789 -v F5ColZB 12.95 -v F5ColXC 11.6397 -v F5ColYC 5.789 -v F5ColZC 12.95 -v F5ColZG  15.95
# -T surface object TopMod 2  -TMod energy 1 "5E-9"
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' pipe3011.x


flat: ess
	pstudy-nowrap -i xml/flat-master.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r  -lowMod None -topModCooling Onion -topPreCooling Onion -x flat.xml -T point free  489.016592713104 0 14.45 -TMod energy 5 "5E-9 2E-8 1E-7 4E-7  1.0 3E+3" -TMod time 5 "1E+8" --validCheck 1000 flat
#-TW -WP 0 0 14.45 flat
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' flat1.x


# www.vtk.org
flatvtk: ess
	pstudy-nowrap -i xml/flat-master.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r -lowMod None -topModCooling Onion -topPreCooling Onion -x flat.xml  -vtk -M -MA -100 -100 -100 -MB 100 100 100  -MN 200 200 200   flatvtk
# -a   axis reorientation
# -M  vtk on the material selection
#-TW -WP 0 0 14.45 flat
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' flat1.x


# flatpipe: ess
# 	pstudy-nowrap -i xml/flat-master.xml -setup
# 	/bin/mv -f case001/inp flat.xml  && rm -fr case001
# 	@rm -f flat1.x
# 	./ess -r -lowMod None -topModCooling Onion -topPreCooling None -x flat.xml  -v ProtonBeamWindowActive 0 -v TopBeRefVoidThick 0 -v TopPreHeight3 0 -v TopPreThick3 0  -v WheelMainShaftMat 26000 -v WheelCladShaftMat 26000  -T point free  489.016592713104 0 14.45 -TMod energy 5 "5E-9 2E-8 1E-7 4E-7  1.0 3E+3" -TMod time 5 "1E+8"  --validCheck 1000 flatpipe
# 	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' flatpipe1.x
# 	perl -pi -e "s/nps 20000/stop f5 1E-2/" flatpipe1.x

flatpipe: ess
	pstudy-nowrap -i xml/flat-master.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r -lowMod None -topModCooling Onion -topPreCooling Onion -x flat.xml -T point free  489.016592713104 0 14.45 -TMod energy 5 "5E-9 2E-8 1E-7 4E-7  1.0 3E+3" -TMod time 5 "1E+8"  --validCheck 1000 flatpipe
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' flatpipe1.x
	perl -pi -e "s/nps 20000/stop f5 1E-2/" flatpipe1.x

flatpipe10cm: ess
	pstudy-nowrap -i xml/flat-master10cm.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r -lowMod None -topModCooling Onion -topPreCooling Onion -x flat.xml -T point free  489.016592713104 0 17.95 -TMod energy 5 "5E-9 2E-8 1E-7 4E-7  1.0 3E+3" -TMod time 5 "1E+8"  --validCheck 1000 flatpipe10cm
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' flatpipe10cm1.x
	perl -pi -e "s/nps 20000/stop f5 1E-2/" flatpipe10cm1.x

flatpipe10cmFOM: ess
	pstudy-nowrap -i xml/flat-master10cm.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r -lowMod None -topModCooling Onion1 -topPreCooling Onion1 -x flat.xml -T point free  489.016592713104 0 17.95 -TMod energy 5 "5E-9 2E-8 1E-7 4E-7  1.0 3E+3" -TMod time 5 "1E+8"  --validCheck 1000 flatpipe10cm
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' flatpipe10cm1.x
	perl -pi -e "s/nps 20000/stop f5 1E-2/" flatpipe10cm1.x

flatpipeTube: ess
	pstudy-nowrap -i xml/flat-master.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r -lowMod Tube -v LowBeRefHeight 40 -v LowBeRefWallMat 13000 -v LowBeRefRefMat 4010 -v TopBeRefRefMat 4010 -v LowBeRefVoidThick 1 -v LowAFlightHeight 3 -v LowAFlightWidth 3 -v LowAFlightXStep 0 -v LowAFlightAngleXY1 1  -v LowAFlightAngleXY2 1 -topModCooling Onion -topPreCooling Onion -v BulkMat1 26010 -v BulkMat3 26010 -x flat.xml -T point free  -178 0 -17.8 -TMod energy 5 "5E-9 2E-8 1E-7 4E-7  1.0 3E+3" -TMod time 5 "1E+8"  --validCheck 1000 flatpipeTube
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' flatpipeTube1.x
	perl -pi -e "s/nps 20000/stop f5 1E-2/" flatpipeTube1.x

thermalGroove: ess
	pstudy-nowrap -i xml/flat-master.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r -lowMod GroovedTube  -v TubeModZStep 3.7 -v TubeModTubeHeight 6 -v TubeModTubeDepth 6 -v GrooveNGrooves 1 -v GrooveLength 10 -v GrooveHeight 5.9 -v GrooveWidth 0.3 -v GrooveDY 1 -v GrooveZStep 3 -v TubeModNLayers 2 -v TubeModPreWingThick 15.3 -v TubeModPreWallThick 0 -v TubeModModMat 1011 -v TubeModMat1 1011 -v TubeModDepth 0 -v TubeModHeight 0 -v TubeModPreHeight1 0 -v TubeModPreHeight2 0 -v TubeModPreDepth 0 -v LowBeRefYStep 0 -v LowBeRefHeight 50 -v LowBeRefWallMat 13000 -v LowBeRefRefMat 4010 -v TopBeRefRefMat 4010 -v LowBeRefVoidThick 1 -v LowAFlightHeight 3 -v LowAFlightWidth 3 -v LowAFlightXStep 0 -v LowAFlightAngleXY1 1  -v LowAFlightAngleXY2 1 -topModCooling Onion -topPreCooling Onion -v BulkMat1 26010 -v BulkMat3 26010 -x flat.xml  --validCheck 1000 thermalGroove
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' thermalGroove1.x
#	perl -pi -e "s/nps 20000/stop f5 1E-2/" thermalGroove1.x

thermal1: ess
	pstudy-nowrap -i xml/flat-master.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r -lowMod Tube  -v TubeModZStep 3.7 -v TubeModTubeHeight 6 -v TubeModTubeDepth 6 -v GrooveNGrooves 1 -v GrooveLength 10 -v GrooveHeight 5.9 -v GrooveWidth 0.3 -v GrooveDY 1 -v GrooveZStep 3 -v TubeModNLayers 2 -v TubeModPreWingThick 15.3 -v TubeModPreWallThick 0 -v TubeModModMat 1011 -v TubeModMat1 1011 -v TubeModDepth 0 -v TubeModHeight 0 -v TubeModPreHeight1 0 -v TubeModPreHeight2 0 -v TubeModPreDepth 0 -v LowBeRefYStep 0 -v LowBeRefHeight 50 -v LowBeRefWallMat 13000 -v LowBeRefRefMat 4010 -v TopBeRefRefMat 4010 -v LowBeRefVoidThick 1 -v LowAFlightHeight 3 -v LowAFlightWidth 3 -v LowAFlightXStep 0 -v LowAFlightAngleXY1 1  -v LowAFlightAngleXY2 1 -topModCooling Onion -topPreCooling Onion -v BulkMat1 26010 -v BulkMat3 26010 -x flat.xml  --validCheck 1000 thermal
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' thermal1.x

thermal: ess
	pstudy-nowrap -i xml/flat-master.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r -lowMod Thermal -v LowBeRefYStep 0 -v LowBeRefHeight 50 -v LowBeRefWallMat 13000 -v LowBeRefRefMat 4010 -v TopBeRefRefMat 4010 -v LowBeRefVoidThick 1 -v LowAFlightHeight 3 -v LowAFlightWidth 3 -v LowAFlightXStep 0 -v LowAFlightAngleXY1 1  -v LowAFlightAngleXY2 1 -topModCooling Onion -topPreCooling Onion -v BulkMat1 26010 -v BulkMat3 26010 -x flat.xml  --validCheck 1000 thermal
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' thermal1.x

thermalcyl: ess
	pstudy-nowrap -i xml/flat-master-water-cylinder.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r -waterCylinder On -thermalCylMod Off -topModCooling Onion -topPreCooling Onion -x flat.xml thermalcyl
#              -v F5ColLength 840 -v F5ColXStep 717 -v F5ColYStep 451 -v F5ColZStep -15.41 -v F5ColXB 10.58 -v F5ColYB 10.37 -v F5ColZB -20.57 -v F5ColXC 14.31 -v F5ColYC 3.88 -v F5ColZC -20.57 -v F5ColZG  -10.25 thermalcyl
#              -v F5ColLength 50 -v F5ColXStep 100 -v F5ColYStep 0 -v F5ColZStep 0 -v F5ColXB 0 -v F5ColYB 1 -v F5ColZB -1 -v F5ColXC 0 -v F5ColYC -1 -v F5ColZC -1 -v F5ColZG 1 \
              --validCheck 1000 thermalcyl
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' thermalcyl1.x

m3a: ess
	pstudy-nowrap -i xml/m3a.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r -topMod M3A -topWaterDisc On -thermalCylMod Off -topModCooling Onion -topPreCooling None  -x flat.xml m3a
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' m3a1.x

butterfly: ess
	pstudy-nowrap -i xml/butterfly.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
#	./ess -r -lowMod Butterfly -topMod Butterfly -nF5 6  -topModCooling Onion  -waterCylinder On  -topPreCooling None -x flat.xml -T surface object Bulk 3 -TMod energy 1 "5E-9 20E-9"  butterfly
	./ess -r -lowMod Butterfly -topMod Butterfly -nF5 6  -topModCooling None  -lowWaterDisc On -topWaterDisc On  -topPreCooling None -x flat.xml -T surface object Bulk 3 -TMod energy 1 "5E-9 20E-9"  butterfly
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' butterfly1.x
	sed -i -e "s/fq1   f d u s m e t c/*c1:n 90 5 0/" butterfly1.x
	remove-weights butterfly1.x

butterflyPIPE: ess
	pstudy-nowrap -i xml/butterfly.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r -nF5 6 -topMod M3A -topModCooling Onion -lowMod Butterfly -waterCylinder Off  -topPreCooling None -x flat.xml  butterfly
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' butterfly1.x
	remove-weights butterfly1.x

butterflyVTK: ess
	pstudy-nowrap -i xml/butterfly.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r -nF5 0 -topMod M3A -topModCooling Onion -lowMod Butterfly -waterCylinder Off  -topPreCooling None -v LowBeRefRefMat 0 -v WaterCylinderRefMat 0 -v LowButterflyFlightLineWallMat 0 -v LowButterflyTopPreType 0 -v LowButterflyPreMat 4010  -x flat.xml -vtk -M -MA -15 -20 -15 -MB 15 20 -7  -MN 500 500 160 butterfly.vtk && du -hsc butterfly.vtk && /bin/mv -f butterfly.vtk ~/Downloads/butterfly.vtk

butterflyVTK1: ess
	pstudy-nowrap -i xml/butterfly.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r -nF5 0 -topMod M3A -topModCooling Onion -lowMod Butterfly -waterCylinder Off  -topPreCooling None -v LowBeRefRefMat 0 -v WaterCylinderRefMat 0 -v LowButterflyFlightLineWallMat 0 -v LowButterflyTopPreType 0 -v LowButterflyPreMat 4010  -x flat.xml -vtk -M -MA -35 -35 -15 -MB 35 35 -7  -MN 500 500 160 butterfly.vtk && du -hsc butterfly.vtk && /bin/mv -f butterfly.vtk ~/Downloads/butterfly.vtk

# an input to test many collimators - inherited from m3a
f5: ess
	pstudy-nowrap -i xml/f5.xml -setup
	/bin/mv -f case001/inp flat.xml  && rm -fr case001
	@rm -f flat1.x
	./ess -r -topMod M3A -waterCylinder Off -thermalCylMod Off -topModCooling None -topPreCooling None  -x flat.xml f5
	sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' f51.x



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
