#include <iostream>
#include <mpi.h>
#include <stdexcept>
#include <string.h>
#include "mdi.h"

//
//This demo example runs a ground-state DFT on a single material system.
//
//Periodic BCC AlNi 2x2x2 16 atom super cell using DFT-FE with ONCV pseudopotentials and Gamma point.
//
//Functionality to update COORDS and CELL is also tested
//

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
  int world_rank;
  MPI_Comm_rank(world_comm, &world_rank);

  // Confirm that the code is being run as a driver
  int role;
  MDI_Get_role(&role);
  if ( role != MDI_DRIVER ) {
    throw std::runtime_error("Must run driver_cxx as a DRIVER");
  }

  // Connect to the engine
  MDI_Comm comm;
  MDI_Accept_communicator(&comm);

  // Confirm that the engine has the @DEFAULT node
  int exists = 1;
  MDI_Check_node_exists("@DEFAULT", comm, &exists);
  if ( exists != 1 ) {
    throw std::runtime_error("The engine does not have the @DEFAULT node.");
  }

  MDI_Check_command_exists("@DEFAULT", "EXIT", comm, &exists);
  if ( exists != 1 ) {
    throw std::runtime_error("The engine does not support the EXIT command.");
  }

  // Determine the name of the engine
  char* engine_name = new char[MDI_NAME_LENGTH];
  MDI_Send_command("<NAME", comm);
  MDI_Recv(engine_name, MDI_NAME_LENGTH, MDI_CHAR, comm);

  if ( world_rank == 0 ) {
    std::cout << " Engine name: " << engine_name << std::endl;
  }


  //Send natoms
  int natoms=16;
  MDI_Send_command(">NATOMS", comm);
  MDI_Send(&natoms, 1, MDI_INT, comm);

  //Send elements
  int elements[natoms]={13,28,13,28,13,28,13,28,13,28,13,28,13,28,13,28};
  MDI_Send_command(">ELEMENTS", comm);
  MDI_Send(&elements, natoms, MDI_INT, comm);

  //Send cell
  double cell[9];
  cell[0]=10.9146877116;cell[1]=0.0;cell[2]=0.0;
  cell[3]=0.0;cell[4]=10.9146877116;cell[5]=0.0;
  cell[6]=0.0;cell[7]=0.0;cell[8]=10.9146877116;  
  MDI_Send_command(">CELL", comm);
  MDI_Send(cell, 9, MDI_DOUBLE, comm);

  //Send dimensions
  int dimensions[3];
  dimensions[0]=2;//periodic
  dimensions[1]=2;//periodic
  dimensions[2]=2;//periodic
  MDI_Send_command(">DIMENSIONS", comm);
  MDI_Send(dimensions, 3, MDI_INT, comm);

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

  MDI_Send_command(">COORDS", comm);
  MDI_Send(coords, 3*natoms, MDI_DOUBLE, comm);



  //compute energy
  double energy;
  MDI_Send_command("<ENERGY", comm);
  MDI_Recv(&energy, 1, MDI_DOUBLE, comm);  
  std::cout<<"DFT free energy: "<<energy<<std::endl;

  //get forces
  double forces[3*natoms];
  MDI_Send_command("<FORCES", comm);
  MDI_Recv(forces, 3*natoms, MDI_DOUBLE, comm); 
  for (int i=0; i<natoms; i++)
     std::cout<<"Atomid: "<<i<<", x: "<<forces[3*i+0]<<", y: "<<forces[3*i+1]<<", z: "<<forces[3*i+2]<<std::endl;

  //update coordinates
  coords[0]=0.1;
  coords[1]=0.1;
  coords[2]=0.1;
  MDI_Send_command(">COORDS", comm);
  MDI_Send(coords, 3*natoms, MDI_DOUBLE, comm);

  //compute energy
  MDI_Send_command("<ENERGY", comm);
  MDI_Recv(&energy, 1, MDI_DOUBLE, comm);  
  std::cout<<"DFT free energy after COORDS update: "<<energy<<std::endl;

  //update cell
  cell[0]=10.0;cell[1]=0.1;cell[2]=0.0;
  cell[3]=0.0;cell[4]=11.0;cell[5]=0.0;
  cell[6]=0.2;cell[7]=0.5;cell[8]=10.9146877116;  
  MDI_Send_command(">CELL", comm);
  MDI_Send(cell, 9, MDI_DOUBLE, comm);

  //compute energy
  MDI_Send_command("<ENERGY", comm);
  MDI_Recv(&energy, 1, MDI_DOUBLE, comm);  
  std::cout<<"DFT free energy after CELL update: "<<energy<<std::endl;


  delete[] engine_name;
  // Confirm that the engine supports the EXIT command
  // Send the "EXIT" command to the engine
  MDI_Send_command("EXIT", comm);

  // Synchronize all MPI ranks
  MPI_Barrier(world_comm);
  MPI_Finalize();

  return 0;
}
