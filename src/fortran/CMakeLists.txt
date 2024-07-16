add_executable(${ExampleName} ${ExampleName}.F90)

target_link_libraries(${ExampleName} PRIVATE OpenCMISS::libOpenCMISS_Fortran)

install(TARGETS ${ExampleName} DESTINATION bin)

configure_file(constant_rate.xml constant_rate.xml)

configure_file(visualise.cmgui visualise.cmgui)
