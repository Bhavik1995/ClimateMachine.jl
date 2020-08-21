# # Eddy Diffusivity - Mass Flux module

# This tutorial implements the Eddy-Diffusivity Mass-Flux (EDMF)
# parameterization for sub grid scale (SGS) turbulance and convection.
# the model equations and rational are based on:
# Tan et al. (2018); Cohen et al. (2020); Lopez-Gomez et al. (2020)
# and Jia et al. (2020).
# The key predictents of the model are the SGS vertical fluxes of heat moisture
# and momentum, and the cloud fraction in the host model grid.
# The EDMF parameterization divides the host model's grid into $N \ge 2$ subdomain
# that represents a isotropic turbulent enviroment and coherent updraft(s)
# and or downdafts. The parameterization solves prognostic euqations
# for first and second moments in the subdomains and provides both the
# SGS vertical flux and cloud fraction in the host model grid.

# -------------- Subdomains:
# The model has a grid mean (the dycore stats vector),i=1:N updrafts and 
# a single environment subdomain (subscript "0")
# The grid mean is prognostic in first moment and diagnostic in second moment.
# The updrafts are prognostic in first moment and set to zero in second moment.
# The environment is diagnostic in first moment and prognostic in second moment.

# ## Equations solved:
# -------------- First Moment Equations:
#                grid mean, ⟨⋅⟩
# ```math
# \begin{align}
# \frac{∂ ρ}{∂ t} = - ∇ ⋅ (ρu)
# \frac{∂ ρ ⟨u⟩}{∂ t} = - ∇ ⋅ (-ρaK ∇u0 - ρ*MF_{u}) - ∇ ⋅ (ρ⟨u⟩ ⟨u⟩' ) + S_{surface Friction}
# \frac{∂ ρ ⟨e⟩_{int}}{∂ t} = - ∇ ⋅ (-ρaK ∇e_{int,0} - ρ*MF_{e_int}) - ∇ ⋅ (u ρ ⟨e⟩_{int} ) + S_{microphysics}
# \frac{∂ ρ ⟨q⟩_{tot}}{∂ t} = - ∇ ⋅ (-ρaK ∇E_{tot,0} - ρ*MF_{q_tot}) - ∇ ⋅ (u ρ ⟨q⟩_{tot} ) + S_{microphysics}
# \end{align}
# ```
# here the Mass-Flux term for a variable \phi is defined as:
# ```math
# \begin{align}
# MF_ϕ = \sum_{i \ge 1}^{N}(a_i * (\bar{w}_i-\bar{w}0)(\bar{ϕ}_i-\bar{ϕ}0))
# \end{align}
# ```

# and
# K is the Eddy_Diffusivity, given as a function of environmental variables
#
#                i'th updraft equations (no second order flux)
# ```math
# \begin{align}
# \frac{∂ ρa_i}{∂ t}           = - ∂/∂z (ρa_i \bar{w}_i) + (E_{i0} - Δ_{i0})
# \frac{∂ ρa_i \bar{w}_i}{∂ t}       = - ∂/∂z (ρa_i \bar{w}_i \bar{w}_i) + (E_{i0}*\bar{w}_0       - Δ_{i0}*u_i) + ρa_i*b - a_i\frac{∂p^†}{∂z}
# \frac{∂ ρa_i \bar{\theta}_{liq,i}}{∂ t} = - ∂/∂z (ρa_i \bar{w}_i *\bar{\theta}_{liq,i}) + (E_{i0}*\bar{\theta}_{liq,0} - Δ_{i0}*\bar{\theta}_{liq,i}) + ρS_{int,i}
# \frac{∂ ρa_i \bar{q}_{tot,i}}{∂ t} = - ∂/∂z (ρa_i \bar{w}_i *\bar{q}_{tot,i}) + (E_{i0}*\bar{q}_{tot,0} - Δ_{i0}*\bar{q}_{tot,i}) + ρS_{tot,i}
# \end{align}
# ```
# here the subdomain and grid mean buoyancies are defined as:
# ```math
# \begin{align}
# b = \frac{\bar{ρ}_i - ρ_{ref}}{⟨ρ⟩}
# ⟨b⟩ = \frac{⟨ρ⟩ - ρ_{ref}}{⟨ρ⟩}
# \end{align}
# ```
#
#                Environment equations first moment
# a_0 = 1-sum_{i \ge 1}^{N}(a_i)
# \bar{w}_0 = (⟨w⟩-sum_{i \ge 1}^{N}a_i*\bar{w}_i)/a_0
# \bar{\theta}_{liq,0} = (⟨\theta⟩_{liq}-sum_{i \ge 1}^{N}(a_i*\bar{\theta}_{liq,i}))/a_0
# \bar{q}_{tot,0} = (⟨q_tot⟩-sum_{i \ge 1}^{N}(a_i*\bar{q}_{tot,i}))/a_0
#
#                environment equations second moment
# ```math
# \begin{align}
# \frac{∂ ρa_0 \overline{ϕ'ψ'}_0}{∂ t} =  - ∇ ⋅ (-ρa_0⋅K⋅∂\overline{ϕ'ψ'}_0/∂z)
# - ∇ ⋅ (u ρa_0⋅\overline{ϕ'ψ'}_0)   + 2ρa_0⋅K(∂\bar{ϕ}/∂z)(∂\bar{ψ}/∂z)
# + (\sum_{i \ge 1}^{N} (E_{i0}*\overline{ϕ'ψ'}_i) - Δ_{i0}*\overline{ϕ'ψ'}_0)
# + ρa_0⋅D_{ϕ'ψ',0} + ρa_0⋅S_{ϕ'ψ',0}
# \end{align}
# ```
# Ideal gas law and subdomain density
# ```math
# \begin{align}
# ρ_i = ⟨p⟩/R_{m,i} * \bar{T}_i
# \end{align}
# ```


# where
#  * `t`           is time
#  * `z`           is height
#  * `ρ`           is the density
#  * `w`           is the vertical velocity
#  * `e_int`       is the internal energy
#  * `\theta_liq`  is the internal energy
#  * `q_tot`       is the total specific humidity
#  * `K`           is the eddy diffusivity
#  * `b`           is the buoyancy
#  * `E_{i0}`      is the entrainment rate from the environment into i
#  * `Δ_{i0}`      is the detrainment rate from i to the environment
#  * `\overline{ϕ'ψ'}_0` is the environmental covariance of ϕ and ψ
#  * `D`           is a covariance dissipation
#  * `S_{ϕ,i}`     is a source of ϕ in the i'th subdomain

# --------------------- Initial Conditions
# Initial conditions are given for all variables in the grid mean, and subdomain variables assume their grid mean values
# ``
# ------- grid mean:
# ρ = hydrostatic reference state - need to compute that state
# ρu = 0
# ρe_int = convert from input profiles
# ρq_tot = convert from input profiles
# ------- updrafts:
# ρa_i = 0.1/N
# ρau = 0
# ρae_int = a*gm.ρe_int
# ρaq_tot = a*gm.ρq_tot
# ------- environment:
# cld_frac = 0.0
# `ϕ'ψ'` = initial covariance profile
# TKE = initial TKE profile
# ``

# --------------------- Boundary Conditions
#           grid mean
# ``
# surface: ρ =
# z_min: ρu = 0
# z_min: ρe_int = 300*cp ; cp=1000
# z_min: ρq_tot = 0.0
# ``

#           i'th updraft
# ``
# z_min: ρ = 1
# z_min: ρu = 0
# z_min: ρe_int = 302*cp ; cp=1000
# z_min: ρq_tot = 0.0
# ``


#### EDMF model kernels

using Printf
using ClimateMachine.Atmos:
    integral_load_auxiliary_state!,
    integral_set_auxiliary_state!,
    atmos_nodal_update_auxiliary_state!

using ClimateMachine.BalanceLaws:
    UpwardIntegrals,
    indefinite_stack_integral!,
    reverse_indefinite_stack_integral!,
    number_states

import ClimateMachine.BalanceLaws:
    update_auxiliary_state!,
    integral_load_auxiliary_state!,
    integral_set_auxiliary_state!

import ClimateMachine.TurbulenceConvection:
    init_aux_turbconv!,
    turbconv_nodal_update_auxiliary_state!,
    turbconv_boundary_state!,
    turbconv_normal_boundary_flux_second_order!

using ClimateMachine.Thermodynamics: air_pressure, air_density


include(joinpath("helper_funcs", "nondimensional_exchange_functions.jl"))
include(joinpath("helper_funcs", "lamb_smooth_minimum.jl"))
include(joinpath("helper_funcs", "subdomain_statistics.jl"))
include(joinpath("helper_funcs", "diagnose_environment.jl"))
include(joinpath("helper_funcs", "edmf_thermo_states.jl"))
include(joinpath("closures", "entr_detr.jl"))
include(joinpath("closures", "pressure.jl"))
include(joinpath("closures", "mixing_length.jl"))
include(joinpath("closures", "turbulence_functions.jl"))
include(joinpath("closures", "surface_functions.jl"))
# include(joinpath("closures", "micro_phys.jl"))


function vars_state(m::NTuple{N, Updraft}, st::UpwardIntegrals, FT) where {N}
    return Tuple{ntuple(i -> vars_state(m[i], st, FT), N)...}
end

function vars_state(::Updraft, ::UpwardIntegrals, FT)
    @vars(H::FT,)
end

function vars_state(m::EDMF, st::UpwardIntegrals, FT)
    @vars(updraft::vars_state(m.updraft, st, FT))
end

function vars_state(m::NTuple{N, Updraft}, st::Auxiliary, FT) where {N}
    return Tuple{ntuple(i -> vars_state(m[i], st, FT), N)...}
end

function vars_state(::Updraft, ::Auxiliary, FT)
    @vars(
        buoyancy::FT,
        updraft_top::FT,
        a::FT,
        ε_dyn::FT,
        δ_dyn::FT,
        ε_trb::FT,
        ε_δ::FT,
        T::FT,
        H::FT,
        H_integ::FT,
        θ_liq::FT,
        q_tot::FT,
        w::FT,
    )
end

function vars_state(::Environment, ::Auxiliary, FT)
    @vars(T::FT, cld_frac::FT, buoyancy::FT)
end

function vars_state(m::EDMF, st::Auxiliary, FT)
    @vars(
        environment::vars_state(m.environment, st, FT),
        updraft::vars_state(m.updraft, st, FT)
    )
end

function vars_state(::Updraft, ::Prognostic, FT)
    @vars(ρa::FT, ρaw::FT, ρaθ_liq::FT, ρaq_tot::FT,)
end

function vars_state(::Environment, ::Prognostic, FT)
    @vars(ρatke::FT, ρaθ_liq_cv::FT, ρaq_tot_cv::FT, ρaθ_liq_q_tot_cv::FT,)
end

function vars_state(m::NTuple{N, Updraft}, st::Prognostic, FT) where {N}
    return Tuple{ntuple(i -> vars_state(m[i], st, FT), N)...}
end


function vars_state(m::EDMF, st::Prognostic, FT)
    @vars(
        environment::vars_state(m.environment, st, FT),
        updraft::vars_state(m.updraft, st, FT)
    )
end

function vars_state(::Updraft, ::Gradient, FT)
    @vars(w::FT,)
end

function vars_state(::Environment, ::Gradient, FT)
    @vars(
        θ_liq::FT,
        q_tot::FT,
        w::FT,
        tke::FT,
        θ_liq_cv::FT,
        q_tot_cv::FT,
        θ_liq_q_tot_cv::FT,
        θv::FT,
        e::FT,
    )
end


function vars_state(m::NTuple{N, Updraft}, st::Gradient, FT) where {N}
    return Tuple{ntuple(i -> vars_state(m[i], st, FT), N)...}
end


function vars_state(m::EDMF, st::Gradient, FT)
    @vars(
        environment::vars_state(m.environment, st, FT),
        updraft::vars_state(m.updraft, st, FT)
    )
end


function vars_state(m::NTuple{N, Updraft}, st::GradientFlux, FT) where {N}
    return Tuple{ntuple(i -> vars_state(m[i], st, FT), N)...}
end

function vars_state(::Updraft, st::GradientFlux, FT)
    @vars(∇w::SVector{3, FT},)
end

function vars_state(::Environment, ::GradientFlux, FT)
    @vars(
        ∇θ_liq::SVector{3, FT},
        ∇q_tot::SVector{3, FT},
        ∇w::SVector{3, FT},
        ∇tke::SVector{3, FT},
        ∇θ_liq_cv::SVector{3, FT},
        ∇q_tot_cv::SVector{3, FT},
        ∇θ_liq_q_tot_cv::SVector{3, FT},
        ∇θv::SVector{3, FT},
        ∇e::SVector{3, FT},
    )

end

function vars_state(m::EDMF, st::GradientFlux, FT)
    @vars(
        S²::FT, # should be conditionally grabbed from atmos.turbulence
        environment::vars_state(m.environment, st, FT),
        updraft::vars_state(m.updraft, st, FT)
    )
end

function init_aux_turbconv!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    aux::Vars,
    geom::LocalGeometry,
) where {FT}

    N_up = n_updrafts(turbconv)

    # Aliases:
    en_a = aux.turbconv.environment
    up_a = aux.turbconv.updraft

    en_a.cld_frac = FT(0)
    en_a.buoyancy = FT(0)

    ntuple(N_up) do i
        up_a[i].buoyancy    = FT(0)
        up_a[i].updraft_top = FT(500)
        up_a[i].θ_liq       = FT(0)
        up_a[i].q_tot       = FT(0)
        up_a[i].w           = FT(0)
    end
end;

# - this method is only called at `t=0`
function init_state_prognostic!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    state::Vars,
    aux::Vars,
    coords,
    t::Real,
) where {FT}

    # # Aliases:
    gm = state
    en = state.turbconv.environment
    up = state.turbconv.updraft
    N_up = n_updrafts(turbconv)
    # GCM setting - Initialize the grid mean profiles of prognostic variables (ρ,e_int,q_tot,u,v,w)
    z = altitude(m, aux)

    # SCM setting - need to have separate cases coded and called from a folder - see what LES does
    # a moist_thermo state is used here to convert the input θ,q_tot to e_int, q_tot profile
    e_int = internal_energy(m, state, aux)

    # Cannot use thermo_state here, since init_aux(::AtmosModel) does not call
    # init_aux(::MoistureModel).
    ts =
        PhaseEquil(m.param_set, e_int, state.ρ, state.moisture.ρq_tot / state.ρ)
    T = air_temperature(ts)
    p = air_pressure(ts)
    q = PhasePartition(ts)
    θ_liq = liquid_ice_pottemp(ts)

    a_up = FT(0.001)
    ntuple(N_up) do i
        up[i].ρa = gm.ρ * a_up
        up[i].ρaw = gm.ρu[3] * a_up
        up[i].ρaθ_liq = gm.ρ * a_up * θ_liq
        up[i].ρaq_tot = gm.moisture.ρq_tot * a_up
    end

    # initialize environment covariance with zero for now
    if (z<=FT(2500))
        en.ρatke = gm.ρ*(FT(1) - z/FT(3000))
    else
        en.ρatke = FT(0)
    end
    en.ρaθ_liq_cv = FT(1e-5) / max(z, FT(10))
    en.ρaq_tot_cv = FT(1e-5) / max(z, FT(10))
    en.ρaθ_liq_q_tot_cv = FT(1e-7) / max(z, FT(10))
    return nothing
end;

# The remaining methods, defined in this section, are called at every
# time-step in the solver by the [`BalanceLaw`](@ref
# ClimateMachine.DGMethods.BalanceLaw) framework.

include("copy_stack_down.jl")

function update_auxiliary_state!(
    dg::DGModel,
    m::AtmosModel,
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
)
    FT = eltype(Q)
    state_auxiliary = dg.state_auxiliary

    if number_states(m, UpwardIntegrals(), FT) > 0
        indefinite_stack_integral!(dg, m, Q, state_auxiliary, t, elems)
        reverse_indefinite_stack_integral!(dg, m, Q, state_auxiliary, t, elems)
    end

    copy_stack_down!(dg, m, m.turbconv, Q, t, elems)

    nodal_update_auxiliary_state!(
        atmos_nodal_update_auxiliary_state!,
        dg,
        m,
        Q,
        t,
        elems,
    )

    # TODO: Remove this hook. This hook was added for implementing
    # the first draft of EDMF, and should be removed so that we can
    # rely on a single vertical element traversal. This hook allows
    # us to compute globally vertical quantities specific to EDMF
    # until we're able to remove them or somehow incorporate them
    # into a higher level hierarchy.
    update_auxiliary_state!(dg, m.turbconv, m, Q, t, elems)

    return true
end

# Compute/update all auxiliary variables at each node. Note that
# - `aux.T` is available here because we've specified `T` in
# `vars_state_auxiliary`
function turbconv_nodal_update_auxiliary_state!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT}

    N_up = n_updrafts(turbconv)
    save_subdomain_temperature!(m, state, aux)

    en_a = aux.turbconv.environment
    up_a = aux.turbconv.updraft
    gm = state
    en = state.turbconv.environment
    up = state.turbconv.updraft

    #  -------------  Compute buoyancies of subdomains
    ts = thermo_state(m, state, aux)
    gm_p = air_pressure(ts)
    ρinv = 1 / gm.ρ
    _grav::FT = grav(m.param_set)

    z = altitude(m, aux)
    ntuple(N_up) do i
        # w_i = state.turbconv.updraft[i].ρaw / state.turbconv.updraft[i].ρa
        ρaw_i = max(state.turbconv.updraft[i].ρaw, 0)
        aux.turbconv.updraft[i].H_integ = max(0, z^10 * ρaw_i)
    end

    ts = thermo_state_en(m, state, aux)
    en_ρ = air_density(ts)
    en_a.buoyancy = -_grav * (en_ρ - aux.ref_state.ρ) * ρinv

    ntuple(N_up) do i
        ts = thermo_state_up(m, state, aux, i)
        ρ_i = air_density(ts)
        up_a[i].buoyancy = -_grav * (ρ_i - aux.ref_state.ρ) * ρinv
        up_a[i].a = up[i].ρa * ρinv
        up_a[i].θ_liq = up[i].ρaθ_liq/up[i].ρa
        up_a[i].q_tot = up[i].ρaq_tot/up[i].ρa
        up_a[i].w = up[i].ρaw/up[i].ρa
    end
    b_gm = grid_mean_b(state, aux, N_up)

    # remove the gm_b from all subdomains
    ntuple(N_up) do i
        up_a[i].buoyancy -= b_gm
    end
    en_a.buoyancy -= b_gm

    ntuple(N_up) do i
        ε_dyn, δ_dyn,ε_trb=
            entr_detr(m, m.turbconv.entr_detr, state, aux, t, i)
        up_a[i].ε_dyn = ε_dyn
        up_a[i].δ_dyn = ε_trb[i]
        up_a[i].ε_trb = δ_dyn[i]
        up_a[i].ε_δ = ε_dyn - δ_dyn
    end

end;

enforce_unit_bounds(x::FT) where {FT} = clamp(x, FT(1e-3), FT(1 - 1e-3))
enforce_positivity(x::FT) where {FT} = max(x, FT(0))

# Since we have second-order fluxes, we must tell `ClimateMachine` to compute
# the gradient of `ρcT`. Here, we specify how `ρcT` is computed. Note that
#  - `transform.ρcT` is available here because we've specified `ρcT` in
#  `vars_state_gradient`
function compute_gradient_argument!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT}
    N_up = n_updrafts(turbconv)
    z = altitude(m, aux)

    # Aliases:
    up_t = transform.turbconv.updraft
    en_t = transform.turbconv.environment
    gm = state
    up = state.turbconv.updraft
    en = state.turbconv.environment

    ntuple(N_up) do i
        up_t[i].w = up[i].ρaw / up[i].ρa
    end
    _grav::FT = grav(m.param_set)
    ts = thermo_state(m, state, aux)

    ρinv = 1 / gm.ρ
    en_area = environment_area(state, aux, N_up)
    en_w = environment_w(state, aux, N_up)
    ts_en = thermo_state_en(m, state, aux)
    en_θ_liq = liquid_ice_pottemp(ts_en)
    en_q_tot = total_specific_humidity(ts_en)

    # populate gradient arguments
    # en_t.θ_liq = en_θ_liq
    en_t.q_tot = en_q_tot
    en_t.w = en_w

    en_t.tke = enforce_positivity(en.ρatke) / (en_area * gm.ρ)
    en_t.θ_liq_cv = enforce_positivity(en.ρaθ_liq_cv) / (en_area * gm.ρ)
    en_t.q_tot_cv = enforce_positivity(en.ρaq_tot_cv) / (en_area * gm.ρ)
    en_t.θ_liq_q_tot_cv = en.ρaθ_liq_q_tot_cv / (en_area * gm.ρ)

    en_t.θv = virtual_pottemp(ts)
    gm_p = air_pressure(ts)
    ts_en = thermo_state_en(m, state, aux)
    e_kin = FT(1 // 2) * ((gm.ρu[1] * ρinv)^2 + (gm.ρu[2] * ρinv)^2 + en_w^2)
    en_t.e = total_energy(e_kin, _grav * z, ts_en)
end;

# Specify where in `diffusive::Vars` to store the computed gradient from
function compute_gradient_flux!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT}
    N_up = n_updrafts(turbconv)

    # # Aliases:
    gm = state
    gm_d = diffusive
    up_d = diffusive.turbconv.updraft
    tc_d = diffusive.turbconv
    up_∇t = ∇transform.turbconv.updraft
    en_d = diffusive.turbconv.environment
    en_∇t = ∇transform.turbconv.environment
    tc_∇t = ∇transform.turbconv

    ntuple(N_up) do i
        up_d[i].∇w = up_∇t[i].w
    end

    ρinv = 1 / gm.ρ
    # negative signs here as we have a '-' sign in BL form leading to + K∂ϕ/∂z on the RHS
    # first moment grid mean coming from environment gradients only
    en_d.∇θ_liq = en_∇t.θ_liq
    en_d.∇q_tot = en_∇t.q_tot
    en_d.∇w = en_∇t.w
    # second moment env cov
    en_d.∇tke = en_∇t.tke
    en_d.∇θ_liq_cv = en_∇t.θ_liq_cv
    en_d.∇q_tot_cv = en_∇t.q_tot_cv
    en_d.∇θ_liq_q_tot_cv = en_∇t.θ_liq_q_tot_cv

    en_d.∇θv = en_∇t.θv
    en_d.∇e = en_∇t.e

    tc_d.S² = ∇transform.u[3, 1]^2 + ∇transform.u[3, 2]^2 + en_d.∇w[3]^2
end;

# We have no sources, nor non-diffusive fluxes.

function turbconv_source!(
    m::AtmosModel{FT},
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
) where {FT}

    turbconv = m.turbconv
    N_up = n_updrafts(turbconv)

    # Aliases:
    gm = state
    en = state.turbconv.environment
    up = state.turbconv.updraft
    gm_s = source
    gm_d = diffusive
    en_s = source.turbconv.environment
    up_s = source.turbconv.updraft
    en_d = diffusive.turbconv.environment
    up_a = aux.turbconv.updraft

    # grid mean sources - I think that large scale subsidence in
    #            doubly periodic domains should be applied here
    ε_trb = MArray{Tuple{N_up}, FT}(zeros(FT, N_up))
    ε_dyn = MArray{Tuple{N_up}, FT}(zeros(FT, N_up))
    δ_dyn = MArray{Tuple{N_up}, FT}(zeros(FT, N_up))

    # # get environment values for e, q_tot , u[3]
    _grav::FT = grav(m.param_set)
    ρinv = 1 / gm.ρ
    en_a = environment_area(state, aux, N_up)
    en_w = environment_w(state, aux, N_up)
    ts_en = thermo_state_en(m, state, aux)
    en_θ_liq = liquid_ice_pottemp(ts_en)
    en_q_tot = total_specific_humidity(ts_en)
    tke_env = enforce_positivity(en.ρatke) * ρinv / en_a
    ts = thermo_state(m, state, aux)
    gm_θ_liq = liquid_ice_pottemp(ts)

    ntuple(N_up) do i
        # upd vars
        w_i = up[i].ρaw / up[i].ρa
        ρa_i = enforce_unit_bounds(up[i].ρa)

        # first moment sources - for now we compute these as aux variable
        ε_dyn[i], δ_dyn[i], ε_trb[i] =
            entr_detr(m, m.turbconv.entr_detr, state, aux, t, i)
        # δ_dyn[i] = FT(0)
        dpdz, dpdz_tke_i = perturbation_pressure(
            m,
            m.turbconv.pressure,
            state,
            diffusive,
            aux,
            t,
            i,
        )

        # entrainment and detrainment
        up_s[i].ρa += up[i].ρaw * (ε_dyn[i] - δ_dyn[i])
        up_s[i].ρaw +=
            up[i].ρaw * (
                (ε_dyn[i] + ε_trb[i]) * en_w -
                (δ_dyn[i] + ε_trb[i]) * up[i].ρaw / ρa_i
            )
        up_s[i].ρaθ_liq +=
            up[i].ρaw * (
                (ε_dyn[i] + ε_trb[i]) * en_θ_liq -
                (δ_dyn[i] + ε_trb[i]) * up[i].ρaθ_liq / ρa_i
            )
        up_s[i].ρaq_tot +=
            up[i].ρaw * (
                (ε_dyn[i] + ε_trb[i]) * en_q_tot -
                (δ_dyn[i] + ε_trb[i]) * up[i].ρaq_tot / ρa_i
            )

        # add buoyancy and perturbation pressure in subdomain w equation
        up_s[i].ρaw += up[i].ρa * (up_a[i].buoyancy - dpdz)
        # microphysics sources should be applied here

        # environment second moments:

        # covariances entrainment sources from the i'th updraft
        # need to compute e_int in updraft and gridmean for entrainment
        # -- if ϕ'ψ' is tke and ϕ,ψ are both w than a factor 0.5 appears in the εt and δ terms
        # Covar_Source      +=  ρaw⋅δ⋅(ϕ_up-ϕ_en) ⋅ (ψ_up-ψ_en)
        #                     + ρaw⋅εt⋅ [(ϕ_up-⟨ϕ⟩)⋅(ψ_up-ψ_en) + (ψ_up-⟨ψ⟩)⋅(ϕ_up-ϕ_en)]
        #                     - ρaw⋅(ε+εt)⋅ϕ'ψ'

        en_s.ρatke += (
            up[i].ρaw *
            δ_dyn[i] *
            (up[i].ρaw / ρa_i - en_w) *
            (up[i].ρaw / ρa_i - en_w) *
            FT(0.5) +
            up[i].ρaw *
            ε_trb[i] *
            (en_w - gm.ρu[3] * ρinv) *
            (en_w - up[i].ρaw / ρa_i) -
            up[i].ρaw * (ε_dyn[i] + ε_trb[i]) * tke_env
        )

        en_s.ρaθ_liq_cv += (
            up[i].ρaw *
            δ_dyn[i] *
            (up[i].ρaθ_liq / ρa_i - en_θ_liq) *
            (up[i].ρaθ_liq / ρa_i - en_θ_liq) +
            up[i].ρaw *
            ε_trb[i] *
            (en_θ_liq - gm_θ_liq) *
            (en_θ_liq - up[i].ρaθ_liq / ρa_i) +
            up[i].ρaw *
            ε_trb[i] *
            (en_θ_liq - gm_θ_liq) *
            (en_θ_liq - up[i].ρaθ_liq / ρa_i) -
            up[i].ρaw * (ε_dyn[i] + ε_trb[i]) * en.ρaθ_liq_cv
        )

        en_s.ρaq_tot_cv += (
            up[i].ρaw *
            δ_dyn[i] *
            (up[i].ρaq_tot / ρa_i - en_q_tot) *
            (up[i].ρaq_tot / ρa_i - en_q_tot) +
            up[i].ρaw *
            ε_trb[i] *
            (en_q_tot - gm.moisture.ρq_tot * ρinv) *
            (en_q_tot - up[i].ρaq_tot / ρa_i) +
            up[i].ρaw *
            ε_trb[i] *
            (en_q_tot - gm.moisture.ρq_tot * ρinv) *
            (en_q_tot - up[i].ρaq_tot / ρa_i) -
            up[i].ρaw * (ε_dyn[i] + ε_trb[i]) * en.ρaq_tot_cv
        )

        en_s.ρaθ_liq_q_tot_cv += (
            up[i].ρaw *
            δ_dyn[i] *
            (up[i].ρaθ_liq / ρa_i - en_θ_liq) *
            (up[i].ρaq_tot / ρa_i - en_q_tot) +
            up[i].ρaw *
            ε_trb[i] *
            (en_θ_liq - gm_θ_liq) *
            (en_q_tot - up[i].ρaq_tot / ρa_i) +
            up[i].ρaw *
            ε_trb[i] *
            (en_q_tot - gm.moisture.ρq_tot * ρinv) *
            (en_θ_liq - up[i].ρaθ_liq / ρa_i) -
            up[i].ρaw * (ε_dyn[i] + ε_trb[i]) * en.ρaθ_liq_q_tot_cv
        )

        # pressure tke source from the i'th updraft
        en_s.ρatke += up[i].ρa * dpdz_tke_i
    end
    l_mix    = mixing_length(m, m.turbconv.mix_len, state, diffusive, aux, t, δ_dyn, ε_trb)
    K_eddy = m.turbconv.mix_len.c_m * l_mix * sqrt(tke_env)
    Shear² = diffusive.turbconv.S²

    # second moment production from mean gradients (+ sign here as we have + S in BL form)
    # production from mean gradient           - Dissipation
    en_s.ρatke +=
        gm.ρ *
        en_a *
        (
            K_eddy * Shear² -
            m.turbconv.mix_len.c_d * sqrt(tke_env) / l_mix * tke_env
        )
    en_s.ρaθ_liq_cv +=
        gm.ρ *
        en_a *
        (
            K_eddy * en_d.∇θ_liq[3] * en_d.∇θ_liq[3] -
            m.turbconv.mix_len.c_d * sqrt(tke_env) / l_mix * en.ρaθ_liq_cv
        )
    en_s.ρaq_tot_cv +=
        gm.ρ *
        en_a *
        (
            K_eddy * en_d.∇q_tot[3] * en_d.∇q_tot[3] -
            m.turbconv.mix_len.c_d * sqrt(tke_env) / l_mix * en.ρaq_tot_cv
        )
    en_s.ρaθ_liq_q_tot_cv +=
        gm.ρ *
        en_a *
        (
            K_eddy * en_d.∇θ_liq[3] * en_d.∇q_tot[3] -
            m.turbconv.mix_len.c_d * sqrt(tke_env) / l_mix *
            en.ρaθ_liq_q_tot_cv
        )
    # covariance microphysics sources should be applied here
end;

# # in the EDMF first order (advective) fluxes exist only in the grid mean (if <w> is nonzero) and the uprdafts
function flux_first_order!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT}

    # # Aliases:
    gm = state
    up = state.turbconv.updraft
    up_f = flux.turbconv.updraft
    N_up = n_updrafts(turbconv)

    # positive sign here as we have a '-' sign in BL form leading to - ∂ρwϕ/∂z on the RHS
    # updrafts
    # in GCM implementation we need to think about grid mean advection
    ρinv = 1 / gm.ρ
    ntuple(N_up) do i
        ρa_i = enforce_unit_bounds(up[i].ρa)
        up_f[i].ρa = up[i].ρaw
        w = up[i].ρaw / ρa_i
        up_f[i].ρaw = up[i].ρaw * w
        up_f[i].ρaθ_liq = w * up[i].ρaθ_liq
        up_f[i].ρaq_tot = w * up[i].ρaq_tot
    end
end;

# in the EDMF second order (diffusive) fluxes exist only in the grid mean and the environment
function flux_second_order!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
) where {FT}
    N_up = n_updrafts(turbconv)

    # Aliases:
    gm = state
    up = state.turbconv.updraft
    en = state.turbconv.environment
    gm_f = flux
    up_f = flux.turbconv.updraft
    en_f = flux.turbconv.environment
    en_d = diffusive.turbconv.environment
    up_a = aux.turbconv.updraft
    ρinv = FT(1) / gm.ρ
    _grav::FT = grav(m.param_set)
    z = altitude(m, aux)

    ε_dyn = MArray{Tuple{N_up}, FT}(zeros(FT, N_up))
    δ_dyn = MArray{Tuple{N_up}, FT}(zeros(FT, N_up))
    ε_trb = MArray{Tuple{N_up}, FT}(zeros(FT, N_up))

    ntuple(N_up) do i
        # for now we compute these as aux variable
        ε_dyn[i], δ_dyn[i], ε_trb[i] = entr_detr(m, m.turbconv.entr_detr, state, aux, t, i)
        ε_dyn[i] = up_a[i].ε_dyn
        δ_dyn[i] = up_a[i].δ_dyn
        ε_trb[i] = up_a[i].ε_trb
        # δ_dyn[i] = FT(0)
    end
    l_mix = mixing_length(m, turbconv.mix_len, state, diffusive, aux, t, δ_dyn, ε_trb)
    en_area = environment_area(state, aux, N_up)
    tke_env = enforce_positivity(en.ρatke) / en_area * ρinv
    K_eddy = m.turbconv.mix_len.c_m * l_mix * sqrt(tke_env)

    # we are adding the massflux term here as it is part of the total flux:
    #total flux(ϕ) =   diffusive_flux(ϕ)  +   massflux(ϕ)
    #   ⟨w ⃰ ϕ ⃰ ⟩   = - a_0 K_eddy⋅∂ϕ/∂z + ∑ a_i(w_i-⟨w⟩)(ϕ_i-⟨ϕ⟩)

    massflux_e = FT(0)
    e_int = internal_energy(m, state, aux)
    ts = thermo_state(m, state, aux)
    gm_p = air_pressure(ts)

    ntuple(N_up) do i
        ts = thermo_state_up(m, state, aux, i)
        e_kin =
            FT(1 // 2) * (
                (gm.ρu[1] * ρinv)^2 +
                (gm.ρu[2] * ρinv)^2 +
                (up[i].ρaw / up[i].ρa)^2
            )
        up_e = total_energy(e_kin, _grav * z, ts)
        ρa_i = enforce_unit_bounds(up[i].ρa)
        massflux_e +=
            up[i].ρa *
            ρinv *
            (gm.ρe * ρinv - up_e) *
            (gm.ρu[3] * ρinv - up[i].ρaw / ρa_i)
    end

    massflux_q_tot = sum(ntuple(
        i ->
            up[i].ρa *
            ρinv *
            (gm.moisture.ρq_tot * ρinv - up[i].ρaq_tot / up[i].ρa) *
            (gm.ρu[3] * ρinv - up[i].ρaw / enforce_unit_bounds(up[i].ρa)),
        N_up,
    ))

    massflux_w = sum(ntuple(
        i ->
            up[i].ρa *
            ρinv *
            (gm.ρu[3] * ρinv - up[i].ρaw / up[i].ρa) *
            (gm.ρu[3] * ρinv - up[i].ρaw / enforce_unit_bounds(up[i].ρa)),
        N_up,
    ))

    # # update grid mean flux_second_order
    ρe_sgs_flux = -gm.ρ * en_area * K_eddy * en_d.∇e[3] + massflux_e
    ρq_tot_sgs_flux = -gm.ρ * en_area * K_eddy * en_d.∇q_tot[3] + massflux_q_tot
    ρu_sgs_flux = -gm.ρ * en_area * K_eddy * en_d.∇w[3] + massflux_w

    # gm_f.ρe              += SVector{3,FT}(0,0,ρe_sgs_flux)
    # gm_f.moisture.ρq_tot += SVector{3,FT}(0,0,ρq_tot_sgs_flux)
    # gm_f.ρu              += SMatrix{3, 3, FT, 9}(
    #     0, 0, 0,
    #     0, 0, 0,
    #     0, 0, ρu_sgs_flux,
    # )

    # env second moment flux_second_order
    en_f.ρatke = -gm.ρ * en_area * K_eddy * en_d.∇tke[3]
    en_f.ρaθ_liq_cv = -gm.ρ * en_area * K_eddy * en_d.∇θ_liq_cv[3]
    en_f.ρaq_tot_cv = -gm.ρ * en_area * K_eddy * en_d.∇q_tot_cv[3]
    en_f.ρaθ_liq_q_tot_cv = -gm.ρ * en_area * K_eddy * en_d.∇θ_liq_q_tot_cv[3]
end;

# ### Boundary conditions

# Second-order terms in our equations, ``∇⋅(G)`` where ``G = α∇ρcT``, are
# internally reformulated to first-order unknowns.
# Boundary conditions must be specified for all unknowns, both first-order and
# second-order unknowns which have been reformulated.

# The boundary conditions for `ρaq_tot` (first order unknown)
function turbconv_boundary_state!(
    nf,
    bc::EDMFBCs,
    m::AtmosModel{FT},
    state⁺::Vars,
    aux⁺::Vars,
    n⁻,
    state⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    state_int::Vars,
    aux_int::Vars,
) where {FT}

    turbconv = m.turbconv
    N_up = n_updrafts(turbconv)
    up = state⁺.turbconv.updraft
    en = state⁺.turbconv.environment
    gm = state⁺
    gm_a = aux⁺
    if bctype == 1 # bottom
        # YAIR - which 'state' should I use here , state⁺ or state⁻  for computation of surface processes
        zLL = altitude(m, aux_int) # or just zLL = FT(20)
        upd_a_surf, upd_θ_liq_surf, upd_q_tot_surf,
        θ_liq_cv, q_tot_cv, θ_liq_q_tot_cv, tke =
            subdomain_surface_values(turbconv.surface,
                                     turbconv,
                                     m,
                                     gm,
                                     gm_a,
                                     zLL)

        ntuple(N_up) do i
            up[i].ρaw = FT(0)
            up[i].ρa = upd_a_surf[i] * gm.ρ
            up[i].ρaθ_liq = up[i].ρa * upd_θ_liq_surf[i]
            up[i].ρaq_tot = up[i].ρa * upd_q_tot_surf[i]
        end
        en_area = environment_area(gm, gm_a, N_up)
        en.ρatke = gm.ρ * en_area * tke
        en.ρaθ_liq_cv = gm.ρ * en_area * θ_liq_cv
        en.ρaq_tot_cv = gm.ρ * en_area * q_tot_cv
        en.ρaθ_liq_q_tot_cv = gm.ρ * en_area * θ_liq_q_tot_cv
    end
end;

# The boundary conditions for second-order unknowns
function turbconv_normal_boundary_flux_second_order!(
    nf,
    bc::EDMFBCs,
    m::AtmosModel{FT},
    fluxᵀn::Vars,
    n⁻,
    state⁻::Vars,
    diff⁻::Vars,
    hyperdiff⁻::Vars,
    aux⁻::Vars,
    state⁺::Vars,
    diff⁺::Vars,
    hyperdiff⁺::Vars,
    aux⁺::Vars,
    bctype,
    t,
    _...,
) where {FT}

    turbconv = m.turbconv
    N = n_updrafts(turbconv)
    up_f = fluxᵀn.turbconv.updraft
    en_f = fluxᵀn.turbconv.environment
    if bctype == 2 # top
        ntuple(N_up) do i
            up_f[i].ρaw      = -n⁻ * FT(0)
            up_f[i].ρa       = -n⁻ * FT(0)
            up_f[i].ρaθ_liq  = -n⁻ * FT(0)
            up_f[i].ρaq_tot  = -n⁻ * FT(0)
        end
        en_f.∇tke            = -n⁻ * FT(0)
        en_f.∇e_int_cv       = -n⁻ * FT(0)
        en_f.∇q_tot_cv       = -n⁻ * FT(0)
        en_f.∇e_int_q_tot_cv = -n⁻ * FT(0)
    end
end;

function integral_load_auxiliary_state!(
    turbconv::EDMF,
    bl::AtmosModel,
    integ::Vars,
    state::Vars,
    aux::Vars,
)
    z = altitude(bl, aux)
    N_up = n_updrafts(turbconv)

    ntuple(N_up) do i
        # w_i = state.turbconv.updraft[i].ρaw / state.turbconv.updraft[i].ρa
        ρaw_i = max(state.turbconv.updraft[i].ρaw, 0)
        integ.turbconv.updraft[i].H = max(0, z^10 * ρaw_i)
        # println("--------------------------------------")
        # @show z^10, ρaw_i
        # if integ.turbconv.updraft[i].H < 0
        #     @show z^10, ρaw_i
        #     @show i, integ.turbconv.updraft[i].H
        #     error("Bad integ.turbconv.updraft[i].H in integral_load_auxiliary_state!")
        # else
        #     println("Good values in integral_load_auxiliary_state!")
        #     @show z^10, ρaw_i
        #     @show i, integ.turbconv.updraft[i].H
        # end
    end
    return nothing
end

function integral_set_auxiliary_state!(
    turbconv::EDMF,
    bl::AtmosModel,
    aux::Vars,
    integ::Vars,
)
    N_up = n_updrafts(turbconv)
    z = altitude(bl, aux)

    ntuple(N_up) do i
        # if integ.turbconv.updraft[i].H < 0
        #     @show z, i, integ.turbconv.updraft[i].H
        #     # error("Bad integ.turbconv.updraft[i].H in integral_set_auxiliary_state!")
        # else
        #     println("Good values in integral_set_auxiliary_state!")
        #     @show z, i, integ.turbconv.updraft[i].H
        # end
        # aux.∫dz.turbconv.updraft[i].H = (integ.turbconv.updraft[i].H)^(1/10)
    end
    return nothing
end
