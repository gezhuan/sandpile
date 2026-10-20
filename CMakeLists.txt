
########################################################################
# MechSys - Open Library for Mechanical Systems
# Copyright (C) 2010 Sergio Galindo
########################################################################

cmake_minimum_required(VERSION 3.25)
project(myproject)

INCLUDE(UseMECHSYS.cmake)

# CUDA architectures supported by Westlake University HPC:
# A100 = 80, A40 = 86, L40S = 89
# Override with -DMECHSYS_CUDA_ARCHITECTURES="80;86;89"
SET(MECHSYS_CUDA_ARCHITECTURES "80;86;89" CACHE STRING
    "CUDA architectures for WL HPC")

SET(PROGS
    # tlbmdem01
    # tlbmdem02
)

SET(CUDAPROGS
    # tlbmdem_cu_01
    # SGCBORUNGPU
    # Sheargpu
    suspensioncu
)

FOREACH(var ${PROGS})
    ADD_EXECUTABLE(${var} "${var}.cpp")
    TARGET_LINK_LIBRARIES(${var} ${LIBS})

    SET_TARGET_PROPERTIES(${var} PROPERTIES
        COMPILE_FLAGS "${FLAGS}"
        LINK_FLAGS "${LFLAGS}"
    )
ENDFOREACH()

FOREACH(var ${CUDAPROGS})
    ADD_EXECUTABLE(${var} "${var}.cu")
    TARGET_LINK_LIBRARIES(${var} ${CUDA_LIBRARIES})

    SET_TARGET_PROPERTIES(${var} PROPERTIES
        COMPILE_FLAGS "${FLAGS} -DUSE_CUDA"
        LINK_FLAGS "${LFLAGS}"
        CUDA_ARCHITECTURES "${MECHSYS_CUDA_ARCHITECTURES}"
    )
ENDFOREACH()
