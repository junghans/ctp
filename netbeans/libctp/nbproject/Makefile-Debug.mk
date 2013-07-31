#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/484457893/segment.o \
	${OBJECTDIR}/_ext/484457893/parallelpaircalc.o \
	${OBJECTDIR}/_ext/238600121/neighborlist.o \
	${OBJECTDIR}/_ext/484457893/qmapplication.o \
	${OBJECTDIR}/_ext/484457893/apolarsite.o \
	${OBJECTDIR}/_ext/238600121/idft.o \
	${OBJECTDIR}/_ext/238600121/tdump.o \
	${OBJECTDIR}/_ext/238600121/izindo.o \
	${OBJECTDIR}/_ext/484457893/qmpair.o \
	${OBJECTDIR}/_ext/484457893/gaussian.o \
	${OBJECTDIR}/_ext/484457893/polarsite.o \
	${OBJECTDIR}/_ext/484457893/qmdatabase.o \
	${OBJECTDIR}/_ext/484457893/topology.o \
	${OBJECTDIR}/_ext/484457893/version_nb.o \
	${OBJECTDIR}/_ext/484457893/qmtopology.o \
	${OBJECTDIR}/_ext/484457893/fragment.o \
	${OBJECTDIR}/_ext/238600121/rates.o \
	${OBJECTDIR}/_ext/238600121/emultipole_stdal.o \
	${OBJECTDIR}/_ext/484457893/orbitals.o \
	${OBJECTDIR}/_ext/484457893/molecule.o \
	${OBJECTDIR}/_ext/238600121/stateserver.o \
	${OBJECTDIR}/_ext/484457893/statesaversqlite.o \
	${OBJECTDIR}/_ext/484457893/qmnblist.o \
	${OBJECTDIR}/_ext/238600121/sandbox.o \
	${OBJECTDIR}/_ext/484457893/calculatorfactory.o \
	${OBJECTDIR}/_ext/238600121/einternal.o \
	${OBJECTDIR}/_ext/238600121/eoutersphere.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a

${OBJECTDIR}/_ext/484457893/segment.o: ../../src/libctp/segment.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/segment.o ../../src/libctp/segment.cc

${OBJECTDIR}/_ext/484457893/parallelpaircalc.o: ../../src/libctp/parallelpaircalc.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/parallelpaircalc.o ../../src/libctp/parallelpaircalc.cc

${OBJECTDIR}/_ext/238600121/neighborlist.o: ../../src/libctp/calculators/neighborlist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/neighborlist.o ../../src/libctp/calculators/neighborlist.cc

${OBJECTDIR}/_ext/484457893/qmapplication.o: ../../src/libctp/qmapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmapplication.o ../../src/libctp/qmapplication.cc

${OBJECTDIR}/_ext/484457893/apolarsite.o: ../../src/libctp/apolarsite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/apolarsite.o ../../src/libctp/apolarsite.cc

${OBJECTDIR}/_ext/238600121/idft.o: ../../src/libctp/calculators/idft.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/idft.o ../../src/libctp/calculators/idft.cc

${OBJECTDIR}/_ext/238600121/tdump.o: ../../src/libctp/calculators/tdump.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/tdump.o ../../src/libctp/calculators/tdump.cc

${OBJECTDIR}/_ext/238600121/izindo.o: ../../src/libctp/calculators/izindo.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/izindo.o ../../src/libctp/calculators/izindo.cc

${OBJECTDIR}/_ext/484457893/qmpair.o: ../../src/libctp/qmpair.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmpair.o ../../src/libctp/qmpair.cc

${OBJECTDIR}/_ext/484457893/gaussian.o: ../../src/libctp/gaussian.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/gaussian.o ../../src/libctp/gaussian.cc

${OBJECTDIR}/_ext/484457893/polarsite.o: ../../src/libctp/polarsite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/polarsite.o ../../src/libctp/polarsite.cc

${OBJECTDIR}/_ext/484457893/qmdatabase.o: ../../src/libctp/qmdatabase.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmdatabase.o ../../src/libctp/qmdatabase.cc

${OBJECTDIR}/_ext/484457893/topology.o: ../../src/libctp/topology.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/topology.o ../../src/libctp/topology.cc

${OBJECTDIR}/_ext/484457893/version_nb.o: ../../src/libctp/version_nb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/version_nb.o ../../src/libctp/version_nb.cc

${OBJECTDIR}/_ext/484457893/qmtopology.o: ../../src/libctp/qmtopology.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmtopology.o ../../src/libctp/qmtopology.cc

${OBJECTDIR}/_ext/484457893/fragment.o: ../../src/libctp/fragment.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/fragment.o ../../src/libctp/fragment.cc

${OBJECTDIR}/_ext/238600121/rates.o: ../../src/libctp/calculators/rates.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/rates.o ../../src/libctp/calculators/rates.cc

${OBJECTDIR}/_ext/238600121/emultipole_stdal.o: ../../src/libctp/calculators/emultipole_stdal.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/emultipole_stdal.o ../../src/libctp/calculators/emultipole_stdal.cc

${OBJECTDIR}/_ext/484457893/orbitals.o: ../../src/libctp/orbitals.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/orbitals.o ../../src/libctp/orbitals.cc

${OBJECTDIR}/_ext/484457893/molecule.o: ../../src/libctp/molecule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/molecule.o ../../src/libctp/molecule.cc

${OBJECTDIR}/_ext/238600121/stateserver.o: ../../src/libctp/calculators/stateserver.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/stateserver.o ../../src/libctp/calculators/stateserver.cc

${OBJECTDIR}/_ext/484457893/statesaversqlite.o: ../../src/libctp/statesaversqlite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/statesaversqlite.o ../../src/libctp/statesaversqlite.cc

${OBJECTDIR}/_ext/484457893/qmnblist.o: ../../src/libctp/qmnblist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmnblist.o ../../src/libctp/qmnblist.cc

${OBJECTDIR}/_ext/238600121/sandbox.o: ../../src/libctp/calculators/sandbox.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/sandbox.o ../../src/libctp/calculators/sandbox.cc

${OBJECTDIR}/_ext/484457893/calculatorfactory.o: ../../src/libctp/calculatorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/calculatorfactory.o ../../src/libctp/calculatorfactory.cc

${OBJECTDIR}/_ext/238600121/einternal.o: ../../src/libctp/calculators/einternal.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/einternal.o ../../src/libctp/calculators/einternal.cc

${OBJECTDIR}/_ext/238600121/eoutersphere.o: ../../src/libctp/calculators/eoutersphere.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.c) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/eoutersphere.o ../../src/libctp/calculators/eoutersphere.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
