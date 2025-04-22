foreach( LIBRARY IN LISTS LIBRARIES_S_L
                          LIBRARIES_P_L )
  if( ${LIBRARY} MATCHES qlapack-.* )

    set_module_path_directory( ${LIBRARY} )

    target_sources( ${LIBRARY}
        PRIVATE
           device_utils.f90
           dgemm_tn.f90
           dgemv.f90
           dgemvn.f90
           dgemvt.f90
           dger.f90
           dger2.f90
           getrfi.f90
           zgemm_tn.f90
           zgemv.f90
           zgemvn.f90
           zgemvt.f90
           zgeru.f90
           zger2.f90
           zdotc.f
    )
    
  string(REGEX REPLACE "^qlapack" "common" COMMON_LIBRARY_TO_LINK ${LIBRARY})

    target_link_libraries( ${LIBRARY}
        PRIVATE
        ${COMMON_LIBRARY_TO_LINK}
        )

    if( ${LIBRARY} MATCHES .*-nogpu-.* )
      set_target_properties( ${LIBRARY} PROPERTIES NOFLAG "_OFFLOAD;_CUSOLVER;_CUBLAS" )
    endif()

  endif()
endforeach()
