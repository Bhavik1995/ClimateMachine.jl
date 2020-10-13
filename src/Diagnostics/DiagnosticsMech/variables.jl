# TODO:
# - scalars and others from #1596, especially core (conditional)
# - diagnostics from a specialized kernel, e.g. vorticity

"""
    States

Composite of the various states, used as a parameter to diagnostic
collection functions.
"""
struct States{PS,GFS,AS,TS}
    prognostic::PS
    gradient_flux::GFS
    auxiliary::AS
    thermodynamic::TS
end

"""
    AbstractIntermediates

Base type for a composite that will be generated.
"""
abstract type AbstractIntermediates end

"""
    DiagnosticVar

The base type and interface for diagnostic variables.
"""
abstract type DiagnosticVar end
function dv_name end
function dv_attrib end

# Generate a standardized type name from the diagnostic variable name.
dv_type_name(name) =
    "DV_" *
    string(supertype(getproperty(@__MODULE__, Symbol(name))))[1:3] *
    "_" *
    name

# Default method for variable attributes.
dv_attrib(
    ::Type{ClimateMachineConfigType},
    ::Type{DiagnosticVar},
) = Dict()

# Generate the type and interface functions for a diagnostic variable.
function generate_dv_interface(
    dvtype,
    config_type,
    name,
    units = "",
    long_name = "",
    standard_name = "",
)
    dvtypname = dv_type_name(name)
    attrib_ex = quote end
    if any(a -> a != "", [units, long_name, standard_name])
        attrib_ex = quote
            $(esc(:dv_attrib))(::Type{$config_type}, ::Type{$dvtypname}) =
                OrderedDict(
                    "units" => $units,
                    "long_name" => $long_name,
                    "standard_name" => $standard_name,
                )
        end
    end
    quote
        struct $(esc(dvtypname)) <: $(esc(dvtype)) end
        $(esc(:dv_name))(::Type{$config_type}, ::Type{$dvtypname}) = $name
        $(attrib_ex)
    end
end

# Helper to generate the implementation function for one or more
# diagnostic variables.
function generate_dv_function(
    dvtype,
    config_type,
    names,
    impl,
)
    dvfun = Symbol("dv_", dvtype)
    dvtypname_args = map(
        n -> :(::Type{$n}),
        map(dv_type_name, names),
    )
    @capture(impl, (args_) -> (body_)) ||
        error("Bad implementation for $(names[1])")
    fun_args = map(
        a -> :($(a[1])::$(a[2])),
        map(splitarg, args.args),
    )
    quote
        function $(esc(dvfun))(
            ::Type{$config_type},
            ::Union{$(dvtypname_args...)},
            $(fun_args...),
        )
            $(unblock(body))
        end
    end
end

# Interface to generate an implementation function for one or more
# diagnostic variables.
macro diagnostic_impl(
    impl,
    dvtype,
    config_type,
    names...,
)
    generate_dv_function(
        dvtype,
        config_type,
        names,
        impl,
    )
end

# Diagnostic variable types and interfaces to create diagnostic variables
# of these types.

"""
    IntermediateValue
"""
abstract type IntermediateValue <: DiagnosticVar end
dv_IntermediateValue(
    ::Type{ClimateMachineConfigType},
    ::Union{Type{IntermediateValue}},
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

macro intermediate_value(config_type, name)
    generate_dv_interface(IntermediateValue, config_type, name)
end

macro intermediate_value(impl, config_type, name)
    generate_dv_interface(IntermediateValue, config_type, name)
    generate_dv_function(
        IntermediateValue,
        config_type,
        [name],
        impl,
    )
end

macro intermediate_values(impl, config_type, names...)
    for name in names
        generate_dv_interface(IntermediateValue, config_type, name)
    end
    generate_dv_function(
        IntermediateValue,
        config_type,
        names,
        impl,
    )
end

"""
    PointwiseDiagnostic

A diagnostic with the same dimensions as the `from_grid` chosen in
the diagnostics group.
"""
abstract type PointwiseDiagnostic <: DiagnosticVar end
dv_PointwiseDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Union{Type{PointwiseDiagnostic}},
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

macro pointwise_diagnostic(config_type, name)
    generate_dv_interface(PointwiseDiagnostic, config_type, name)
end

macro pointwise_diagnostic(
    config_type,
    name,
    units,
    long_name,
    standard_name,
)
    generate_dv_interface(
        PointwiseDiagnostic,
        config_type,
        name,
        units,
        long_name,
        standard_name,
    )
end

macro pointwise_diagnostic(
    impl,
    config_type,
    name,
    units,
    long_name,
    standard_name,
)
    generate_dv_interface(
        PointwiseDiagnostic,
        config_type,
        name,
        units,
        long_name,
        standard_name,
    )
    generate_dv_function(
        PointwiseDiagnostic,
        config_type,
        [name],
        impl,
    )
end

macro pointwise_diagnostic_impl(
    impl,
    config_type,
    names...,
)
    generate_dv_function(
        PointwiseDiagnostic,
        config_type,
        names,
        impl,
    )
end

"""
    HorizontalAverage

A horizontal reduction into a single vertical dimension.
"""
abstract type HorizontalAverage <: DiagnosticVar end
dv_HorizontalAverage(
    ::Type{ClimateMachineConfigType},
    ::Type{HorizontalAverage},
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

macro horizontal_average(config_type, name)
    generate_dv_interface(HorizontalAverage, config_type, name)
end

macro horizontal_average(
    config_type,
    name,
    units,
    long_name,
    standard_name,
)
    generate_dv_interface(
        HorizontalAverage,
        config_type,
        name,
        units,
        long_name,
        standard_name,
    )
end

macro horizontal_average(
    impl,
    config_type,
    name,
    units,
    long_name,
    standard_name,
)
    generate_dv_interface(
        HorizontalAverage,
        config_type,
        name,
        units,
        long_name,
        standard_name,
    )
    generate_dv_function(
        HorizontalAverage,
        config_type,
        [name],
        impl,
    )
end

macro horizontal_average_impl(
    impl,
    config_type,
    names...,
)
    generate_dv_function(
        HorizontalAverage,
        config_type,
        names,
        impl,
    )
end


#= TODO
macro scalar_diagnostic(
)
=#

include("atmos_diagnostic_funs.jl")

include("atmos_thermo_vars.jl")
include("atmos_les_diagnostic_vars.jl")
include("atmos_gcm_diagnostic_vars.jl")