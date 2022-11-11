#include <iostream>
#include <mpi.h>
#include <stdexcept>
#include <string.h>
#include <cstdlib>
#include <iomanip>
#include "mdi.h"


void mpi_error(const char* errormsg) {
  std::cerr << errormsg << std::endl;
  MPI_Abort(MPI_COMM_WORLD, 1);
}


int code_for_plugin_instance(void* mpi_comm_ptr, MDI_Comm mdi_comm, void* class_object) {
  MPI_Comm mpi_comm = *(MPI_Comm*) mpi_comm_ptr;
  int my_rank;
  MPI_Comm_rank(mpi_comm, &my_rank);

  // Determine the name of the engine
  char* engine_name = new char[MDI_NAME_LENGTH];
  if ( MDI_Send_command("<NAME", mdi_comm) != 0 ) {
    mpi_error("MDI_Send_command returned non-zero exit code.");
  }
  if ( MDI_Recv(engine_name, MDI_NAME_LENGTH, MDI_CHAR, mdi_comm) != 0 ) {
    mpi_error("MDI_Recv returned non-zero exit code.");
  }

  //Send natoms
  int natoms=16;
  MDI_Send_command(">NATOMS", mdi_comm);
  MDI_Send(&natoms, 1, MDI_INT, mdi_comm);

  //Send elements
  int elements[natoms]={13,28,13,28,13,28,13,28,13,28,13,28,13,28,13,28};
  MDI_Send_command(">ELEMENTS", mdi_comm);
  MDI_Send(&elements, natoms, MDI_INT, mdi_comm);

  //Send cell
  double cell[9];
  cell[0]=10.9146877116;cell[1]=0.0;cell[2]=0.0;
  cell[3]=0.0;cell[4]=10.9146877116;cell[5]=0.0;
  cell[6]=0.0;cell[7]=0.0;cell[8]=10.9146877116;  
  MDI_Send_command(">CELL", mdi_comm);
  MDI_Send(cell, 9, MDI_DOUBLE, mdi_comm);

  //Send dimensions
  int dimensions[3];
  dimensions[0]=2;//periodic
  dimensions[1]=2;//periodic
  dimensions[2]=2;//periodic
  MDI_Send_command(">DIMENSIONS", mdi_comm);
  MDI_Send(dimensions, 3, MDI_INT, mdi_comm);

  //Send MP grid
  int mpgrid[3];
  mpgrid[0]=1;
  mpgrid[1]=1;
  mpgrid[2]=1;
  MDI_Send_command(">MONKHORST-PACK_NPOINTS", mdi_comm);
  MDI_Send(mpgrid, 3, MDI_INT, mdi_comm);  

  //Send MP shift
  double mpshift[3];
  mpshift[0]=0;
  mpshift[1]=0;
  mpshift[2]=0;
  MDI_Send_command(">MONKHORST-PACK_SHIFT", mdi_comm);
  MDI_Send(mpshift, 3, MDI_DOUBLE, mdi_comm);  

  //Send spin polarization
  int spin=0;
  MDI_Send_command(">SPIN_POLARIZATION", mdi_comm);
  MDI_Send(&spin, 1, MDI_INT, mdi_comm);

  //
  //Send cartesian coordinates, origin at cell corner
  //

  //initialize with fractional and then convert to cartesian
  double coordsfrac[3*natoms]={0.000000000000,0.000000000000,0.000000000000,
                                      0.250000000000,0.250000000000,0.250000000000,
                                      0.000000000000,0.000000000000,0.500000000000,
                                      0.250000000000,0.250000000000,0.750000000000,
				                              0.000000000000,0.500000000000,0.000000000000,
                                      0.250000000000,0.750000000000,0.250000000000,
                                      0.000000000000,0.500000000000,0.500000000000,
                                      0.250000000000,0.750000000000,0.750000000000,
                                      0.500000000000,0.000000000000,0.000000000000,
                                      0.750000000000,0.250000000000,0.250000000000,
                                      0.500000000000,0.000000000000,0.500000000000,
                                      0.750000000000,0.250000000000,0.750000000000,
                                      0.500000000000,0.500000000000,0.000000000000,
                                      0.750000000000,0.750000000000,0.250000000000,
                                      0.500000000000,0.500000000000,0.500000000000,
                                      0.750000000000,0.750000000000,0.750000000000};

  double coords[3*natoms]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};                                      

  for (unsigned int i=0;i<natoms;++i)
    for (unsigned int j=0;j<3;++j)      
       for (unsigned int l=0;l<3;++l)
          coords[3*i+j]+=cell[3*l+j]*coordsfrac[3*i+l];

  //for (int i=0; i<natoms; i++)
  //   std::cout<<"coords: "<<i<<", x: "<<coords[3*i+0]<<", y: "<<coords[3*i+1]<<", z: "<<coords[3*i+2]<<std::endl;

  MDI_Send_command(">COORDS", mdi_comm);
  MDI_Send(coords, 3*natoms, MDI_DOUBLE, mdi_comm);



  //compute energy
  double energy;
  MDI_Send_command("<ENERGY", mdi_comm);
  MDI_Recv(&energy, 1, MDI_DOUBLE, mdi_comm);  
  std::cout<<"DFT free energy: "<<std::setprecision(10)<<energy<<std::endl;

  //get forces
  double forces[3*natoms];
  MDI_Send_command("<FORCES", mdi_comm);
  MDI_Recv(forces, 3*natoms, MDI_DOUBLE, mdi_comm); 
  for (int i=0; i<natoms; i++)
     std::cout<<"Force, Atomid: "<<i<<", x: "<<forces[3*i+0]<<", y: "<<forces[3*i+1]<<", z: "<<forces[3*i+2]<<std::endl;

  //update coordinates
  coords[0]=0.1;
  coords[1]=0.1;
  coords[2]=0.1;
  MDI_Send_command(">COORDS", mdi_comm);
  MDI_Send(coords, 3*natoms, MDI_DOUBLE, mdi_comm);

  //compute energy
  MDI_Send_command("<ENERGY", mdi_comm);
  MDI_Recv(&energy, 1, MDI_DOUBLE, mdi_comm);  
  std::cout<<"DFT free energy after COORDS update: "<<std::setprecision(10)<<energy<<std::endl;

  //update cell
  cell[0]=10.0;cell[1]=0.1;cell[2]=0.0;
  cell[3]=0.0;cell[4]=11.0;cell[5]=0.0;
  cell[6]=0.2;cell[7]=0.5;cell[8]=10.9146877116;  
  MDI_Send_command(">CELL", mdi_comm);
  MDI_Send(cell, 9, MDI_DOUBLE, mdi_comm);

  //compute energy
  MDI_Send_command("<ENERGY", mdi_comm);
  MDI_Recv(&energy, 1, MDI_DOUBLE, mdi_comm);  
  std::cout<<"DFT free energy after CELL update: "<<std::setprecision(10)<<energy<<std::endl;

  if ( my_rank == 0 ) {
    std::cout << " Engine name: " << engine_name << std::endl;
  }
  delete[] engine_name;

  // Send the "EXIT" command to the engine
  if ( MDI_Send_command("EXIT", mdi_comm) != 0 ) {
    mpi_error("MDI_Send_command returned non-zero exit code.");
  }

  return 0;
}


int main(int argc, char **argv) {
  int ret;

  // Initialize the MPI environment
  MPI_Comm world_comm;
  MPI_Init(&argc, &argv);

  // Initialize MDI
  ret = MDI_Init(&argc, &argv);
  if ( ret != 0 ) {
    throw std::runtime_error("The MDI library was not initialized correctly.");
  }

  // Check for an "-earlyreturn" option
  bool earlyreturn = false;
  for ( int iarg=1; iarg < argc; iarg++ ) {
    if ( strcmp(argv[iarg],"-earlyreturn") == 0 ) {
      earlyreturn = true;
    }
  }
  
  // If the earlyreturn flag was set, return now
  if ( earlyreturn ) {
    MPI_Finalize();
    return 0;
  }

  // Confirm that MDI was initialized successfully
  int initialized_mdi;
  ret = MDI_Initialized(&initialized_mdi);
  if ( ret != 0 ) {
    throw std::runtime_error("MDI_Initialized failed.");
  }
  if ( ! initialized_mdi ) {
    throw std::runtime_error("MDI not initialized: did you provide the -mdi option?.");
  }

  // Get the correct MPI intra-communicator for this code
  ret = MDI_MPI_get_world_comm(&world_comm);
  if ( ret != 0 ) {
    throw std::runtime_error("MDI_MPI_get_world_comm failed.");
  }

  // Number of ranks that will run the driver
  // This is the number of ranks that will NOT run plugin instances
  // The value of this variable is read from the command-line options
  int driver_nranks = -1;

  // Number of ranks running EACH plugin instance
  // The value of this variable is read from the command-line options
  int plugin_nranks = -1;

  // Name of the plugin to use
  // The value of this variable is read from the command-line options
  char* plugin_name = NULL;

  // Read through all the command line options
  int iarg = 1;
  while ( iarg < argc ) {

    if ( strcmp(argv[iarg],"-driver_nranks") == 0 ) {

      // Ensure that the argument to the -driver_nranks option was provided
      if ( argc-iarg < 2 ) {
        mpi_error("The -driver_nranks argument was not provided.");
      }

      // Set driver_nranks
      char* strtol_ptr;
      driver_nranks = strtol( argv[iarg+1], &strtol_ptr, 10 );
      iarg += 2;

    }
    else if ( strcmp(argv[iarg],"-plugin_nranks") == 0 ) {

      // Ensure that the argument to the -plugin_nranks option was provided
      if ( argc-iarg < 2 ) {
        mpi_error("The -plugin_nranks argument was not provided.");
      }

      // Set driver_nranks
      char* strtol_ptr;
      plugin_nranks = strtol( argv[iarg+1], &strtol_ptr, 10 );
      iarg += 2;

    }
    else if ( strcmp(argv[iarg],"-plugin_name") == 0 ) {

      // Ensure that the argument to the -plugin_name option was provided
      if ( argc-iarg < 2 ) {
        mpi_error("The -plugin_name argument was not provided.");
      }

      // Set driver_nranks
      plugin_name = argv[iarg+1];
      iarg += 2;

    }
    else {
      mpi_error("Unrecognized option.");
    }

  }

  // Verify the value of driver_nranks
  if ( driver_nranks < 0 ) {
    mpi_error("Invalid value for driver_nranks [0, inf).");
  }

  // Verify the value of plugin_nranks
  if ( plugin_nranks <= 0 ) {
    mpi_error("Invalid value for plugin_nranks (0, inf).");
  }

  // Verify that the value of driver_nranks and plugin_nranks is consistent with world_size
  int world_size;
  MPI_Comm_size(world_comm, &world_size);
  if ( (world_size - driver_nranks) % plugin_nranks != 0 ) {
    mpi_error("Invalid values for driver_nranks and plugin_nranks: world_size - driver_nranks must be divisible by plugin_nranks.");
  }

  // Verify the value of plugin_name
  if ( plugin_name == NULL ) {
    mpi_error("Plugin name was not provided.");
  }
  
  // Split world_comm into MPI intra-comms for the driver and each plugin
  MPI_Comm intra_comm;
  int my_rank, color, intra_rank;
  MPI_Comm_rank(world_comm, &my_rank);
  if ( my_rank < driver_nranks ) {
    color = 0;
  }
  else {
    color = ( ( my_rank - driver_nranks ) / plugin_nranks ) + 1;
  }
  MPI_Comm_split(world_comm, color, my_rank, &intra_comm);
  MPI_Comm_rank(intra_comm, &intra_rank);

  if ( color == 0 ) { // Driver intra-comm

    if (intra_rank == 0 ) {
      std::cout << "I am the driver" << std::endl;
    }

  }
  else { // Engine instance intra-comm

    if ( intra_rank == 0 ) {
      std::cout << "I am engine instance: " << color << std::endl;
    }

    ////////////////////////////////////////////////////////
/*
    // Initialize and run an instance of the engine library
    MDI_Comm mdi_comm;
    if ( MDI_Open_plugin(plugin_name,
			   "-mdi \"-name MM -role ENGINE -method LINK\"",
			   &intra_comm,
			   &mdi_comm) != 0 ) {
      mpi_error("MDI_Launch_plugin returned non-zero exit code.");
    }
    char* engine_name = new char[MDI_NAME_LENGTH];
    if ( MDI_Send_command("<NAME", mdi_comm) != 0 ) {
      mpi_error("MDI_Send_command returned non-zero exit code.");
    }
    if ( MDI_Recv(engine_name, MDI_NAME_LENGTH, MDI_CHAR, mdi_comm) != 0 ) {
      mpi_error("MDI_Recv returned non-zero exit code.");
    }
    if ( intra_rank == 0 ) {
      std::cout << "TEST VALUE: " << engine_name << std::endl;
    }
    MDI_Close_plugin(mdi_comm);
*/
    ////////////////////////////////////////////////////////

    // Initialize and run an instance of the engine library
    if ( MDI_Launch_plugin(plugin_name,
			   "-mdi \"-name QM -role ENGINE -method LINK\"",
			   &intra_comm,
			   code_for_plugin_instance,
			   NULL) != 0 ) {
      mpi_error("MDI_Launch_plugin returned non-zero exit code.");
    }
  }

  // Synchronize all MPI ranks
  MPI_Barrier(world_comm);
  MPI_Finalize();

  return 0;
}
