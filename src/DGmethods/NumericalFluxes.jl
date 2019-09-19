module NumericalFluxes

export Rusanov, CentralGradPenalty, CentralNumericalFluxDiffusive

using StaticArrays
import ..DGmethods: BalanceLaw, Grad, Vars, vars_state, vars_diffusive,
                    vars_aux, vars_gradient, boundary_state!, wavespeed,
                    flux_nondiffusive!, flux_diffusive!, diffusive!, num_state,
                    num_gradient, gradvariables!

"""
    GradNumericalPenalty

Any `P <: GradNumericalPenalty` should define methods for:

   diffusive_penalty!(gnf::P, bl::BalanceLaw, σ, n⁻, Q⁻, H⁻, α⁻, Q⁺,
                      H⁺, α⁺, t)
   diffusive_boundary_penalty!(gnf::P, bl::BalanceLaw, l_δ, n⁻, l_H⁻, l_Q⁻,
                               l_α⁻, l_H⁺, l_Q⁺, l_α⁺, bctype, t)

"""
abstract type GradNumericalPenalty end

function diffusive_penalty! end
function diffusive_boundary_penalty! end

"""
    CentralGradPenalty <: GradNumericalPenalty

"""
struct CentralGradPenalty <: GradNumericalPenalty end

function diffusive_penalty!(::CentralGradPenalty, bl::BalanceLaw,
                            σ, n⁻, H⁻, Q⁻, α⁻, H⁺, Q⁺, α⁺, t)
  DFloat = eltype(Q⁻)

  @inbounds begin
    Nᵈ = 3
    ngradstate = num_gradient(bl,DFloat)
    δ = similar(VF, Size(Nᵈ, ngradstate))
    for j = 1:ngradstate, i = 1:Nᵈ
      δ[i, j] = n⁻[i] * (H⁺[j] - H⁻[j]) / 2
    end
    diffusive!(bl,
               Vars{vars_diffusive(bl,DFloat)}(σ),
               Grad{vars_gradient(bl,DFloat)}(δ),
               Vars{vars_state(bl,DFloat)}(Q⁻),
               Vars{vars_aux(bl,DFloat)}(α⁻),
               t)
  end
end

function diffusive_boundary_penalty!(nf::CentralGradPenalty, bl::BalanceLaw,
                                     σ, n⁻, H⁻, Q⁻, α⁻, H⁺, Q⁺, α⁺,
                                     bctype, t, Q1, α1)
  DFloat = eltype(H⁺)

  boundary_state!(nf, bl,
                  Vars{vars_state(bl,DFloat)}(Q⁺),
                  Vars{vars_aux(bl,DFloat)}(α⁺),
                  n⁻,
                  Vars{vars_state(bl,DFloat)}(Q⁻),
                  Vars{vars_aux(bl,DFloat)}(α⁻),
                  bctype, t,
                  Vars{vars_state(bl,DFloat)}(Q1),
                  Vars{vars_aux(bl,DFloat)}(α1))

  gradvariables!(bl,
                 Vars{vars_gradient(bl,DFloat)}(H⁺),
                 Vars{vars_state(bl,DFloat)}(Q⁺),
                 Vars{vars_aux(bl,DFloat)}(α⁺),
                 t)

  diffusive_penalty!(nf, bl, σ, n⁻, H⁻, Q⁻, α⁻, H⁺, Q⁺, α⁺, t)
end


"""
    NumericalFluxNonDiffusive

Any `N <: NumericalFluxNonDiffusive` should define the a method for

    numerical_flux_nondiffusive!(nf::N, bl::BalanceLaw, F, n⁻, Q⁻, α⁻, Q⁺,
                                 α⁺, t)

where
- `F` is the numerical flux array
- `n⁻` is the unit normal
- `Q⁻`/`Q⁺` are the minus/positive state arrays
- `t` is the time

An optional method can also be defined for

    numerical_boundary_flux_nondiffusive!(nf::N, bl::BalanceLaw, F, n⁻, Q⁻,
                                          α⁻, Q⁺, α⁺, bctype, t)

"""
abstract type NumericalFluxNonDiffusive end

function numerical_flux_nondiffusive! end

function numerical_boundary_flux_nondiffusive!(nf::NumericalFluxNonDiffusive,
                                               bl::BalanceLaw,
                                               F::MArray{Tuple{nstate}},
                                               n⁻, Q⁻, α⁻, Q⁺, α⁺, bctype, t,
                                               Q1, α1) where {nstate}
  DFloat = eltype(F)

  boundary_state!(nf, bl,
                  Vars{vars_state(bl,DFloat)}(Q⁺),
                  Vars{vars_aux(bl,DFloat)}(α⁺),
                  n⁻,
                  Vars{vars_state(bl,DFloat)}(Q⁻),
                  Vars{vars_aux(bl,DFloat)}(α⁻),
                  bctype, t,
                  Vars{vars_state(bl,DFloat)}(Q1),
                  Vars{vars_aux(bl,DFloat)}(α1))

  numerical_flux_nondiffusive!(nf, bl, F, n⁻, Q⁻, α⁻, Q⁺, α⁺, t)
end



"""
    Rusanov <: NumericalFluxNonDiffusive

The Rusanov (aka local Lax-Friedrichs) numerical flux.

# Usage

    Rusanov()

Requires a `flux_nondiffusive!` and `wavespeed` method for the balance law.
"""
struct Rusanov <: NumericalFluxNonDiffusive end


function numerical_flux_nondiffusive!(::Rusanov, bl::BalanceLaw, F::MArray,
                                      n⁻, Q⁻, α⁻, Q⁺, α⁺, t)
  DFloat = eltype(F)
  nstate = num_state(bl,DFloat)

  λ⁻ = wavespeed(bl, n⁻,
                 Vars{vars_state(bl,DFloat)}(Q⁻),
                 Vars{vars_aux(bl,DFloat)}(α⁻),
                 t)

  F⁻ = similar(F, Size(3, nstate))
  fill!(F⁻, -zero(eltype(F⁻)))

  flux_nondiffusive!(bl, Grad{vars_state(bl,DFloat)}(F⁻),
                     Vars{vars_state(bl,DFloat)}(Q⁻),
                     Vars{vars_aux(bl,DFloat)}(α⁻), t)

  λ⁺ = wavespeed(bl, n⁻,
                 Vars{vars_state(bl,DFloat)}(Q⁺),
                 Vars{vars_aux(bl,DFloat)}(α⁺),
                 t)

  F⁺ = similar(F, Size(3, nstate))
  fill!(F⁺, -zero(eltype(F⁺)))

  flux_nondiffusive!(bl,
                     Grad{vars_state(bl,DFloat)}(F⁺),
                     Vars{vars_state(bl,DFloat)}(Q⁺),
                     Vars{vars_aux(bl,DFloat)}(α⁺),
                     t)

  λ  =  max(λ⁻, λ⁺)

  @inbounds for s = 1:nstate
    F[s] += 0.5 * (n⁻[1] * (F⁻[1, s] + F⁺[1, s]) +
                   n⁻[2] * (F⁻[2, s] + F⁺[2, s]) +
                   n⁻[3] * (F⁻[3, s] + F⁺[3, s]) +
                   λ * (Q⁻[s] - Q⁺[s]))
  end
end

"""
    NumericalFluxDiffusive

Any `N <: NumericalFluxDiffusive` should define the a method for

    numerical_flux_diffusive!(nf::N, bl::BalanceLaw, F, n⁻, Q⁻, σ⁻, α⁻, Q⁺,
                              σ⁺, α⁺, t)

where
- `F` is the numerical flux array
- `n⁻` is the unit normal
- `Q⁻`/`Q⁺` are the minus/positive state arrays
- `σ⁻`/`σ⁺` are the minus/positive diffusive state arrays
- `α⁻`/`α⁺` are the minus/positive auxiliary state arrays
- `t` is the time

An optional method can also be defined for

    numerical_boundary_flux_diffusive!(nf::N, bl::BalanceLaw, F, n⁻, Q⁻, σ⁻,
                                       α⁻, Q⁺, σ⁺, α⁺, bctype, t)

"""
abstract type NumericalFluxDiffusive end

function numerical_flux_diffusive! end

function numerical_boundary_flux_diffusive!(nf::NumericalFluxDiffusive,
                                            bl::BalanceLaw,
                                            F::MArray{Tuple{nstate}},
                                            n⁻, Q⁻, σ⁻, α⁻, Q⁺, σ⁺, α⁺,
                                            bctype, t, Q1, σ1,
                                            α1) where {nstate}
  DFloat = eltype(F)

  boundary_state!(nf, bl, Vars{vars_state(bl,DFloat)}(Q⁺),
                  Vars{vars_diffusive(bl,DFloat)}(σ⁺),
                  Vars{vars_aux(bl,DFloat)}(α⁺),
                  n⁻,
                  Vars{vars_state(bl,DFloat)}(Q⁻),
                  Vars{vars_diffusive(bl,DFloat)}(σ⁻),
                  Vars{vars_aux(bl,DFloat)}(α⁻),
                  bctype, t,
                  Vars{vars_state(bl,DFloat)}(Q1),
                  Vars{vars_diffusive(bl,DFloat)}(σ1),
                  Vars{vars_aux(bl,DFloat)}(α1))

  numerical_flux_diffusive!(nf, bl, F, n⁻, Q⁻, σ⁻, α⁻, Q⁺, σ⁺, α⁺, t)
end

"""
    CentralNumericalFluxDiffusive <: NumericalFluxDiffusive

The central numerical flux for diffusive terms

# Usage

    CentralNumericalFluxDiffusive()

Requires a `flux_diffusive!` for the balance law.
"""
struct CentralNumericalFluxDiffusive <: NumericalFluxDiffusive end


function numerical_flux_diffusive!(::CentralNumericalFluxDiffusive,
                                   bl::BalanceLaw, F::MArray, n⁻,
                                   Q⁻, σ⁻, α⁻, Q⁺, σ⁺, α⁺, t)
  DFloat = eltype(F)
  nstate = num_state(bl,DFloat)

  F⁻ = similar(F, Size(3, nstate))
  fill!(F⁻, -zero(eltype(F⁻)))

  flux_diffusive!(bl,
                  Grad{vars_state(bl,DFloat)}(F⁻),
                  Vars{vars_state(bl,DFloat)}(Q⁻),
                  Vars{vars_diffusive(bl,DFloat)}(σ⁻),
                  Vars{vars_aux(bl,DFloat)}(α⁻),
                  t)

  F⁺ = similar(F, Size(3, nstate))
  fill!(F⁺, -zero(eltype(F⁺)))

  flux_diffusive!(bl,
                  Grad{vars_state(bl,DFloat)}(F⁺),
                  Vars{vars_state(bl,DFloat)}(Q⁺),
                  Vars{vars_diffusive(bl,DFloat)}(σ⁺),
                  Vars{vars_aux(bl,DFloat)}(α⁺),
                  t)

  @inbounds for s = 1:nstate
    F[s] += 0.5 * (n⁻[1] * (F⁻[1, s] + F⁺[1, s]) +
                   n⁻[2] * (F⁻[2, s] + F⁺[2, s]) +
                   n⁻[3] * (F⁻[3, s] + F⁺[3, s]))
  end
end


end
