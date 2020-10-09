"""
    DiagnosticsGroupParams

Base type for any extra paramters.
"""
abstract type DiagnosticsGroupParams end

"""
    DiagnosticsGroup

Holds a set of diagnostics that share a collection interval, a filename
prefix, an output writer, an interpolation, and any extra parameters.
"""
mutable struct DiagnosticsGroup{DGP <: Union{Nothing, DiagnosticsGroupParams}}
    name::String
    init::Function
    collect::Function
    fini::Function
    interval::String
    out_prefix::String
    writer::AbstractWriter
    interpol::Union{Nothing, InterpolationTopology}
    onetime::Bool
    params::DGP

    DiagnosticsGroup(
        name,
        init,
        collect,
        fini,
        interval,
        out_prefix,
        writer,
        interpol,
        onetime,
        params = nothing,
    ) = new{typeof(params)}(
        name,
        init,
        collect,
        fini,
        interval,
        out_prefix,
        writer,
        interpol,
        onetime,
        params,
    )
end

# `GenericCallbacks` implementations for `DiagnosticsGroup`s
function GenericCallbacks.init!(dgngrp::DiagnosticsGroup, solver, Q, param, t)
    @info @sprintf(
        """
    Diagnostics: %s
        initializing at %8.2f""",
        dgngrp.name,
        t,
    )
    dgngrp.init(dgngrp, t)
    dgngrp.collect(dgngrp, t)
    return nothing
end
function GenericCallbacks.call!(dgngrp::DiagnosticsGroup, solver, Q, param, t)
    @tic diagnostics
    @info @sprintf(
        """
    Diagnostics: %s
        collecting at %8.2f""",
        dgngrp.name,
        t,
    )
    dgngrp.collect(dgngrp, t)
    @toc diagnostics
    return nothing
end
function GenericCallbacks.fini!(dgngrp::DiagnosticsGroup, solver, Q, param, t)
    @info @sprintf(
        """
    Diagnostics: %s
        finishing at %8.2f""",
        dgngrp.name,
        t,
    )
    dgngrp.collect(dgngrp, t)
    dgngrp.fini(dgngrp, t)
    return nothing
end

"""
    OnGridKind

The variables in a diagnostic group are computed from either the DG
grid or an interpolated grid.
"""
abstract type OnGridKind end
struct GridDG <: OnGridKind end
struct GridInterpolated <: OnGridKind end

"""
    @diagnostics_group

Generate the functions needed to establish and use a `DiagnosticsGroup`
containing the specified `DiagnosticVar`s.

# TODO: OnGridKind
"""
macro diagnostics_group(
    name,
    config_type,
    params_type,
    vars...,
)
    setup = generate_setup(name, config_type, params_type)
    init = generate_init(name, vars)
    collect = generate_collect(name, vars)
    fini = generate_fini(name, vars)

    return Expr(
        :block,
        esc(setup),
        esc(init),
        esc(collect),
        esc(fini),
    )
end

include("group_gen.jl")

# Pre-defined `DiagnosticsGroup`s
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
