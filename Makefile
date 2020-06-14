FC=gfortran
switch=-O3 -cpp

SRC_DIR=src
BUILD_DIR=build
BIN_DIR=bin
objects = common_var.o main.o init.o interpolate.o optimization.o Langevin.o tools.o friction_daldrop.o

# standard version
optle: builddirs $(objects)
	$(FC) -o $(BIN_DIR)/optle $(switch) $(addprefix $(BUILD_DIR)/, $(objects))
common_var.mod: common_var.o $(SRC_DIR)/common_var.f90
	$(FC) -o $(BUILD_DIR)/common_var.mod -c $(switch) $(SRC_DIR)/common_var.f90
common_var.o: $(SRC_DIR)/common_var.f90
	$(FC) -o $(BUILD_DIR)/common_var.o -c $(switch) $(SRC_DIR)/common_var.f90
main.o: common_var.mod $(SRC_DIR)/main.f90
	$(FC) -o $(BUILD_DIR)/main.o -c $(switch) $(SRC_DIR)/main.f90
init.o: common_var.mod $(SRC_DIR)/init.f90
	$(FC) -o $(BUILD_DIR)/init.o -c $(switch) $(SRC_DIR)/init.f90
interpolate.o: common_var.mod $(SRC_DIR)/interpolate.f90
	$(FC) -o $(BUILD_DIR)/interpolate.o -c $(switch) $(SRC_DIR)/interpolate.f90
optimization.o: common_var.mod $(SRC_DIR)/optimization.f90
	$(FC) -o $(BUILD_DIR)/optimization.o -c $(switch) $(SRC_DIR)/optimization.f90
Langevin.o: common_var.mod $(SRC_DIR)/Langevin.f90
	$(FC) -o $(BUILD_DIR)/langevin.o -c $(switch) $(SRC_DIR)/Langevin.f90
tools.o: common_var.mod $(SRC_DIR)/tools.f90
	$(FC) -o $(BUILD_DIR)/tools.o -c $(switch) $(SRC_DIR)/tools.f90
friction_daldrop.o: common_var.mod $(SRC_DIR)/friction_daldrop.f90
	$(FC) -o $(BUILD_DIR)/friction_daldrop.o -c $(switch) $(SRC_DIR)/friction_daldrop.f90
builddirs:
	mkdir -p $(BUILD_DIR) $(BIN_DIR)
	cp resources/input $(BIN_DIR)/input
	cp resources/colvar $(BIN_DIR)/colvar
	cp resources/RESTART $(BIN_DIR)/RESTART

# version printing detailed trajectories (use opt_niter 1)
debug: $(objects)
	$(FC) -o $(BIN_DIR)/optle-debug $(switch) -DDEBUG $(addprefix $(BUILD_DIR)/, $(objects))

# clean the mess
clean:
	rm common_var.mod
	rm $(BUILD_DIR)/common_var.mod
	rm $(addprefix $(BUILD_DIR)/, $(objects))