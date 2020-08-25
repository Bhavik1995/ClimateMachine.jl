#### Entrainment-Detrainment kernels
include(joinpath("..", "helper_funcs", "diagnose_environment.jl"))

function entr_detr(
    m::AtmosModel{FT},
    entr::EntrainmentDetrainment,
    state::Vars,
    aux::Vars,
    t::Real,
    i::Int,
) where {FT}

    # Alias convention:
    gm = state
    en = state.turbconv.environment
    up = state.turbconv.updraft
    gm_a = aux
    en_a = aux.turbconv.environment
    up_a = aux.turbconv.updraft

    N_up = n_updrafts(m.turbconv)
    ρinv = 1 / gm.ρ
    up_area = up[i].ρa / gm.ρ
    z = altitude(m, aux)

    sqrt_ϵ = FT(0.0001) #sqrt(eps(FT))
    w_min = FT(0.1)
    # precompute vars
    a_en = environment_area(state, aux, N_up)
    w_en = environment_w(state, aux, N_up)
    w_up = up[i].ρaw / up[i].ρa
    sqrt_tke = sqrt(max(en.ρatke, 0) * ρinv / a_en)
    Δw = w_up - w_en
    if Δw<FT(0)
        Δw = min(Δw, -w_min)
    else
        Δw = max(Δw, w_min)
    end
    if w_up<FT(0)
        w_up = min(w_up, -w_min)
    else
        w_up = max(w_up, w_min)
    end
    Δb = up_a[i].buoyancy - en_a.buoyancy

    D_ε, D_δ, M_δ, M_ε =
        nondimensional_exchange_functions(m, entr, state, aux, t, i)

    # I am commenting this out for now, to make sure there is no slowdown here
    Λ_w = abs(Δb/Δw)
    Λ_tke = entr.c_λ * abs(Δb / (max(en.ρatke*ρinv,0) + w_min))
    λ = lamb_smooth_minimum(SVector(Λ_w, Λ_tke))
    # λ = abs(Δb / Δw)

    # compute limiters
    εt_lim = εt_limiter(w_up, sqrt_ϵ)
    ε_lim = ε_limiter(up_area, sqrt_ϵ)
    δ_lim = δ_limiter(up_area, sqrt_ϵ)
    # compute entrainment/detrainmnet components
    # ε_trb = 2 * up_area * entr.c_t * sqrt_tke /
           # max( (w_up * up_area * up_a[i].updraft_top),FT(1e-4))
    ε_trb =
        2 * up_area * entr.c_t * sqrt_tke /
        max((w_up * up_area * FT(500)), sqrt_ϵ)#* εt_lim
    ε_dyn = λ / w_up * (D_ε + M_ε)*ε_lim
    δ_dyn = λ / w_up * (D_δ + M_δ)*δ_lim

    ε_dyn = min(max(ε_dyn, FT(0)), FT(1))
    δ_dyn = min(max(δ_dyn, FT(0)), FT(1))
    ε_trb = min(max(ε_trb, FT(0)), FT(1))

    return ε_dyn, δ_dyn, ε_trb, D_ε, M_ε, D_δ, M_δ, ε_lim, δ_lim, εt_lim,  λ, Δw
end;

ε_limiter(a_up::FT, ϵ::FT) where {FT} =
    FT(1) + FT(10)*exp(-a_up^2/(2*ϵ)) - exp(-(FT(1)-a_up)^2/(2*ϵ))
δ_limiter(a_up::FT, ϵ::FT) where {FT} =
    FT(1) - exp(-a_up^2/(2*ϵ)) + FT(10)*exp(-(FT(1)-a_up)^2/(2*ϵ))
εt_limiter(w_up::FT, ϵ::FT) where {FT} =
    FT(1) + FT(10)*exp(-w_up^2/(2*ϵ))

# ε_limiter(a_up::FT, ϵ::FT) where {FT} = 1
# δ_limiter(a_up::FT, ϵ::FT) where {FT} = 1
