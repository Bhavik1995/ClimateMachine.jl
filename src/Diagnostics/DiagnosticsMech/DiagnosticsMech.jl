"""
    DiagnosticsMech

This module provides the infrastructure to extract diagnostics from a
ClimateMachine simulation. Two key abstractions are defined: diagnostic
variables and diagnostic groups. The `StdDiagnostics` module makes use of
these to define many standard variables and groups which may be used
directly by experiments. `DiagnosticsMech` may be used by experiments to
define specialized variables and groups.
"""
module DiagnosticsMech

export DiagnosticVar,
    @intermediate_value,
    @intermediate_values,
    @pointwise_diagnostic,
    @pointwise_diagnostic_impl,
    @horizontal_average,
    @horizontal_average_impl,
    DiagnosticsGroup,
    @diagnostics_group,
    DiagnosticsGroupParams

using CUDA
using Dates
using FileIO
using JLD2
using KernelAbstractions
using MacroTools
using MPI
using OrderedCollections
using Printf
using StaticArrays

using ..BalanceLaws
using ..ConfigTypes
using ..DGMethods
using ..GenericCallbacks
using ..Mesh.Interpolation
using ..MPIStateArrays
using ..Spectra
using ..TicToc
using ..VariableTemplates
using ..Writers

using CLIMAParameters
using CLIMAParameters.Planet: planet_radius

# Container to store simulation information necessary for all
# diagnostics groups.
Base.@kwdef mutable struct Diagnostic_Settings
    mpicomm::MPI.Comm = MPI.COMM_WORLD
    param_set::Union{Nothing, AbstractParameterSet} = nothing
    dg::Union{Nothing, DGModel} = nothing
    Q::Union{Nothing, MPIStateArray} = nothing
    starttime::Union{Nothing, String} = nothing
    output_dir::Union{Nothing, String} = nothing
end
const Settings = Diagnostic_Settings()

"""
    init(mpicomm, param_set, dg, Q, starttime, output_dir)

Save the parameters into `Settings`, a container for simulation
information necessary for all diagnostics groups.
"""
function init(
    mpicomm::MPI.Comm,
    param_set::AbstractParameterSet,
    dg::DGModel,
    Q::MPIStateArray,
    starttime::String,
    output_dir::String,
)
    Settings.mpicomm = mpicomm
    Settings.param_set = param_set
    Settings.dg = dg
    Settings.Q = Q
    Settings.starttime = starttime
    Settings.output_dir = output_dir
end

include("helpers.jl")
include("variables.jl")
include("groups.jl")

end # module DiagnosticsMech
