
DATAMODEL=/HOME/ntag/NTag_ToolFramework/DataModel
DATAMODEL=/home/sk_user/ToolAnalysis/DataModel

all: main

main: main.cpp $(DATAMODEL)/EventTrueCaptures.cpp $(DATAMODEL)/Cluster.cpp $(DATAMODEL)/Cluster.h $(DATAMODEL)/EventTrueCaptures.h $(DATAMODEL)/TrueCapture.cpp $(DATAMODEL)/TrueCapture.h $(DATAMODEL)/Particle.cpp $(DATAMODEL)/Particle.h
	g++ -ggdb3 -std=c++11 -lgfortran main.cpp $(DATAMODEL)/EventTrueCaptures.cpp  $(DATAMODEL)/Cluster.cpp $(DATAMODEL)/TrueCapture.cpp $(DATAMODEL)/Particle.cpp $(SKOFL_ROOT)/src/sklib/skprint.o `root-config --cflags --libs` -I $(DATAMODEL) -L $(DATAMODEL)/../lib -lRootDict -o $@
