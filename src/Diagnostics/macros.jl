"""
    DiagnosticVar

A diagnostic variable is an `n`-dimensional field computed from
the various `ClimateMachine` state variables. Diagnostic variables
may be assembled into `DiagnosticsGroup`s.

Standard names are from:
http://cfconventions.org/Data/cf-standard-names/71/build/cf-standard-name-table.html
"""

"""
    DiagnosticVarKind

Kinds of diagnostic variables.
"""
abstract type DiagnosticVarKind end
struct Diagnostic <: DiagnosticVarKind end
struct Variance <: DiagnosticVarKind end
struct Covariance <: DiagnosticVarKind end
struct HorizontalAverage <: DiagnosticVarKind end

"""
    VarLayoutKind

The layout of a diagnostic variable -- may be on the DG grid, or on a
Cartesian grid.
"""
abstract type VarLayoutKind end
struct LayoutDG <: VarLayoutKind end
struct Layout3D <: VarLayoutKind
    x::Int
    y::Int
    z::Int
end
struct Layout2D <: VarLayoutKind
    x::Int
    y::Int
end
struct Layout1D <: VarLayoutKind
    x::Int
end

"""
    DiagnosticVar

A diagnostic variable has a name, attributes (NetCDF), an implementation,
a kind, and a layout. Diagnostic variables may be assembled into groups.
"""
struct DiagnosticVar{F}
    kind::DiagnosticVarKind
    layout::VarLayoutKind
    name::String
    attrib::OrderedDict
    impl::F

    DiagnosticVar(
        impl,
        kind::DiagnosticVarKind,
        layout::VarLayoutKind
        name::String,
        attrib::OrderedDict = OrderedDict(),
    ) = new(kind, layout, name, attrib, impl)
end
const AllDiagnosticVars = OrderedDict{Symbol, DiagnosticVar}()

include("macros.jl")
"""
    @horizontal_average name attrib layout expr
"""
macro horizontal_average(name, attrib, impl)
    quote
        AllDiagnosticVars[Symbol(name)] =
            DiagnosticVar(HorizontalAverage(), name, $attrib, $impl)
    end
end

@horizontal_average(
    u,
    var_attrib("", "", ""),
    () -> conservative.ρu[1],
)
@horizontal_average(
    w_ht_sgs,
    var_attrib("", "", ""),
    (atmos, conservative, gradient_flux, auxiliary, curr_time) -> begin
        ν, D_t, _ = turbulence_tensors(
            atmos,
            conservative,
            gradient_flux,
            auxiliary,
            curr_time,
        )
        d_h_tot = -D_t .* gradient_flux.∇h_tot
        d_h_tot[end]
    end,
)

@intermediate(
    "thermo",
    (atmos, conservative, auxiliary) -> begin
        thermo_state(atmos, conservative, auxiliary)
    end,
)
@diagnostic3d(
    temp,
    var_attrib("", "", ""),
    (atmos, conservative, auxiliary) -> begin
    end,
)

function vars_atmos_les_default_simple(m::AtmosModel, FT)
    @vars begin
        u::FT
        v::FT
        w::FT
        avg_rho::FT             # ρ
        rho::FT                 # ρρ
        temp::FT
        pres::FT
        thd::FT                 # θ_dry
        et::FT                  # e_tot
        ei::FT                  # e_int
        ht::FT
        hi::FT
        w_ht_sgs::FT

        moisture::vars_atmos_les_default_simple(m.moisture, FT)
    end
end
vars_atmos_les_default_simple(::MoistureModel, FT) = @vars()
function vars_atmos_les_default_simple(m::EquilMoist, FT)
    @vars begin
        qt::FT                  # q_tot
        ql::FT                  # q_liq
        qv::FT                  # q_vap
        thv::FT                 # θ_vir
        thl::FT                 # θ_liq
        w_qt_sgs::FT
    end
end

@diagnostic_group "AtmosLESDefault" vars_atmos_les_default_simple(m, FT) 

