"""
    StdDiagnostics

This module defines many standard diagnostic variables and groups that may
be used directly by experiments.
"""
module StdDiagnostics

using ..DiagnosticsMech

export setup_atmos_default_diagnostics,
    setup_atmos_core_diagnostics,
    setup_atmos_default_perturbations,
    setup_atmos_refstate_perturbations,
    setup_atmos_turbulence_stats,
    setup_atmos_mass_energy_loss,
    setup_dump_state_diagnostics,
    setup_dump_aux_diagnostics,
    setup_dump_spectra_diagnostics

# Pre-defined diagnostic variables
include("atmos_les_diagnostic_vars.jl")
include("atmos_gcm_diagnostic_vars.jl")

# Pre-defined diagnostics groups
include("atmos_les_default.jl")
include("atmos_gcm_default.jl")
include("atmos_les_core.jl")
include("atmos_les_default_perturbations.jl")
include("atmos_refstate_perturbations.jl")
include("atmos_turbulence_stats.jl")
include("atmos_mass_energy_loss.jl")
include("dump_init.jl")
include("dump_state.jl")
include("dump_aux.jl")
include("dump_spectra.jl")

end # module StdDiagnostics
