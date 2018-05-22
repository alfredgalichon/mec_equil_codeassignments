using DataFrames
using Roots
using CSV

include("surgeHelper.jl")
include("surgeSolve.jl")
include("surgeTax.jl")

###############
# Model Setup #
###############

# Parameters
precision = 1e-5;

####################
# gradient descent #
####################

GDstats = gradientDescent(CYZ, my, UXZ, nx, precision)

##################
# Newton Descent #
##################

NDstats = newtonMethod(CYZ, my, UXZ, nx, precision)

###############################
# Coordinate Descent (Jacobi) #
###############################

CDstats = coordDescentJacobi(CYZ,my,UXZ,nx, precision)

#####################################
# Coordinate Descent (Gauss-Seidel) #
#####################################

CDGSstats = coordDescentGaussSeidel(CYZ,my,UXZ,nx,precision)

############################
# Modification: New Supply #
############################

CDstatstax = coordDescentJacobitax(CYZ,my,UXZ,nx,precision)
CDGSstatstax, pCDGStax = coordDescentGaussSeideltax(CYZ,my,UXZ,nx,precision)

################
# Save to File #
################

CSV.write("surgeStats.csv",DataFrame(hcat(GDstats, NDstats, CDstats, CDGSstats)))
CSV.write("surgeTaxStats.csv",DataFrame(hcat(CDstatstax, CDGSstatstax)))
