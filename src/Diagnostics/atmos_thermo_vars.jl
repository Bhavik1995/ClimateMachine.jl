using ..Atmos
using ..ConfigTypes
using ..Thermodynamics

@intermediate_values(
    AtmosConfigType,
    "ts",
    "temp",
    "pres",
    "θ_dry",
    "e_int",
    "h_tot",
    "h_int",
) do (
    atmos::AtmosModel,
    states::States,
    curr_time,
    intermediates,
)
    e_tot = states.prognostic.ρe / states.prognostic.ρ
    ts = recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    temp = air_temperature(ts)
    pres = air_pressure(ts)
    θ_dry = dry_pottemp(ts)
    e_int = internal_energy(ts)
    h_tot = total_specific_enthalpy(ts, e_tot)
    h_int = specific_enthalpy(ts)
end

@intermediate_values(
    AtmosConfigType,
    "q_liq",
    "q_ice",
    "q_vap",
    "θ_vir",
    "θ_liq_ice",
    "has_condensate",
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states::States,
    curr_time,
    intermediates,
)
    q_liq = liquid_specific_humidity(intermediates.ts)
    q_ice = ice_specific_humidity(intermediates.ts)
    q_vap = vapor_specific_humidity(intermediates.ts)
    θ_vir = virtual_pottemp(intermediates.ts)
    θ_liq_ice = liquid_ice_pottemp(intermediates.ts)
    has_condensate = has_condensate(intermediates.ts)
end
