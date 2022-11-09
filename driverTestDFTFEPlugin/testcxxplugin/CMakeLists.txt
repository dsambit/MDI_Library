# Check for MPI

if ( NOT ( mpi STREQUAL "OFF") )
   find_package(MPI)
endif()
if( NOT MPI_FOUND )
   if( mpi STREQUAL "ON" )
      message( WARNING "Could not find MPI.  Compiling without MPI support." )
   endif()
   configure_file(${CMAKE_CURRENT_SOURCE_DIR}/STUBS_MPI/mpi.h ${CMAKE_CURRENT_BINARY_DIR}/STUBS_MPI/mpi.h COPYONLY)
   set(mpi "OFF")
endif()



# Macro to link target to MPI

macro(link_against_mpi _TAR)
   if( mpi STREQUAL "ON" )

      #include MPI
      string_to_list(MPI_CXX_COMPILE_OPTIONS   "${MPI_CXX_COMPILE_FLAGS}")
      string_to_list(MPI_CXX_LINK_OPTIONS      "${MPI_CXX_LINK_FLAGS}")

      target_include_directories(${_TAR} PRIVATE ${MPI_CXX_INCLUDE_PATH})
      target_compile_options(${_TAR} PRIVATE ${MPI_CXX_COMPILE_OPTIONS})
      target_link_libraries(${_TAR} PRIVATE ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_OPTIONS})

   elseif( mpi STREQUAL "OFF" )

      message( "Compiling without MPI." )
      target_include_directories(${_TAR} PRIVATE ${CMAKE_BINARY_DIR}/tests/MDI_Test_Codes/STUBS_MPI/)

   else()

      message( FATAL_ERROR "Value of mpi not recognized. Accepted values are: ON; OFF." )

   endif()
endmacro(link_against_mpi _TAR)



# Link everything to MDI

link_directories( ${mdi_lib_path} )
include_directories( ${mdi_include_path} )



if ( use_CXX )

   if ( test_drivers )
      add_subdirectory(driver_plug_cxx)
      #add_subdirectory(driver_repeat_cxx)
      #add_subdirectory(driver_serial_cxx)
      #add_subdirectory(driver_plug_cxx)
   endif()

   #if ( test_engines )
      #add_subdirectory(engine_cxx)
      # Can't compile the i-PI test on Windows
      #if(NOT WIN32)
      #   add_subdirectory(engine_ipi_cxx)
      #endif()
   #endif()

endif()
