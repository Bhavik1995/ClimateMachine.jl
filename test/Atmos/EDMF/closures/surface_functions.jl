#### Surface model kernels
## --- revert to use compute_buoyancy_flux in SurfaceFluxes.jl ---|

using Statistics

function subdomain_surface_values(
    m::SurfaceModel,
    turbconv::EDMF{FT},
    atmos::AtmosModel{FT},
    state::Vars,
    aux::Vars,
    zLL::FT,
) where {FT}

    turbconv = atmos.turbconv
    N_up = n_updrafts(turbconv)
    gm = state
    en = state.turbconv.environment
    up = state.turbconv.updraft
    ts = thermo_state(atmos, state, aux)
    gm_p = air_pressure(ts)
    θ_liq = liquid_ice_pottemp(ts)
    q = PhasePartition(ts)
    _cp_m = cp_m(atmos.param_set, q)
    lv = latent_heat_vapor(ts)
    Π = exner(ts)
    ρinv = 1 / gm.ρ
    surface_scalar_coeff = turbconv.surface.scalar_coeff

    θ_liq_surface_flux = m.surface_shf / Π / _cp_m
    q_tot_surface_flux = m.surface_lhf / lv
    # these value should be given from the SurfaceFluxes.jl once it is merged
    oblength = -FT(100)
    ustar = FT(0.28)

    if oblength < 0
        θ_liq_cv =
            4 * (θ_liq_surface_flux * θ_liq_surface_flux) / (ustar * ustar) *
            (1 - FT(8.3) * zLL / oblength)^(-FT(2) / FT(3))
        q_tot_cv =
            4 * (q_tot_surface_flux * q_tot_surface_flux) / (ustar * ustar) *
            (1 - FT(8.3) * zLL / oblength)^(-FT(2) / FT(3))
        θ_liq_q_tot_cv =
            4 * (θ_liq_surface_flux * q_tot_surface_flux) / (ustar * ustar) *
            (1 - FT(8.3) * zLL / oblength)^(-FT(2) / FT(3))
        tke = ustar * ustar * (FT(3.75) + cbrt(zLL / oblength * zLL / oblength))
    else
        θ_liq_cv =
            4 * (θ_liq_surface_flux * θ_liq_surface_flux) / (ustar * ustar)
        q_tot_cv =
            4 * (q_tot_surface_flux * q_tot_surface_flux) / (ustar * ustar)
        θ_liq_q_tot_cv =
            4 * (θ_liq_surface_flux * q_tot_surface_flux) / (ustar * ustar)
        tke = FT(3.75) * ustar * ustar
    end

    upd_a_surf = MArray{Tuple{N_up}, FT}(zeros(FT, N_up))
    upd_θ_liq_surf = MArray{Tuple{N_up}, FT}(zeros(FT, N_up))
    upd_q_tot_surf = MArray{Tuple{N_up}, FT}(zeros(FT, N_up))
    # ntuple(N_up) do i
    #     upd_a_surf[i] = FT(m.a_surf / N_up)
    #     e_int = internal_energy(atmos, state, aux)
    #     ts = PhaseEquil(
    #         atmos.param_set,
    #         e_int,
    #         state.ρ,
    #         state.moisture.ρq_tot / state.ρ,
    #     )
    #     gm_θ_liq = liquid_ice_pottemp(ts)
    #     upd_θ_liq_surf[i] =
    #         (gm_θ_liq + surface_scalar_coeff[i] * sqrt(max(θ_liq_cv, 0)))
    #     upd_q_tot_surf[i] = (
    #         gm.moisture.ρq_tot * ρinv +
    #         surface_scalar_coeff[i] * sqrt(max(q_tot_cv, 0))
    #     )
    # end
    return upd_a_surf, upd_θ_liq_surf, upd_q_tot_surf, θ_liq_cv, q_tot_cv, θ_liq_q_tot_cv, tke
end;

function percentile_bounds_mean_norm(
    low_percentile::FT,
    high_percentile::FT,
    n_samples::IT,
) where {FT <: Real, IT}
    Random.seed!(15)
    x = rand(Normal(), n_samples)
    xp_low = quantile(Normal(), low_percentile)
    xp_high = quantile(Normal(), high_percentile)
    filter!(y -> xp_low < y < xp_high, x)
    return Statistics.mean(x)
end
