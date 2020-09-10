module Spectra

export power_spectrum, power_spectrum_zonal, power_spectrum_gcm_2d

using MPI
using ClimateMachine
using ..BalanceLaws
using ..ConfigTypes
using ..DGMethods
using ..Mesh.Interpolation
using ..Mesh.Grids
using ..MPIStateArrays
using ..VariableTemplates
using ..Writers

include("power_spectrum.jl")
include("power_spectrum_gcm.jl")

end
