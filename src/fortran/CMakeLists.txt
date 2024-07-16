add_executable(darcy_static_fortran darcy_static.F90)

target_link_libraries(darcy_static_fortran PRIVATE OpenCMISS::libOpenCMISS_Fortran)

install(TARGETS darcy_static_fortran DESTINATION bin)
