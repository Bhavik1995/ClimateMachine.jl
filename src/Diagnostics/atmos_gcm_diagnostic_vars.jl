@pointwise_diagnostic(
    "u",
    AtmosGCMConfigType,
    GridInterpolated,
    "m s^-1",
    "zonal wind",
    "eastward_wind",
) do (atmos::AtmosModel, states, curr_time)
    states.prognostic.ρu[1] / states.prognostic.ρ
end

@pointwise_diagnostic(
    "v",
    AtmosGCMConfigType,
    GridInterpolated,
    "m s^-1",
    "meridional wind",
    "northward_wind",
) do (atmos::AtmosModel, states, curr_time)
    states.prognostic.ρu[2] / states.prognostic.ρ
end

@pointwise_diagnostic(
    "w",
    AtmosGCMConfigType,
    GridInterpolated,
    "m s^-1",
    "vertical wind",
    "upward_air_velocity",
) do (atmos::AtmosModel, states, curr_time)
    states.prognostic.ρu[3] / states.prognostic.ρ
end

@pointwise_diagnostic(
    "rho",
    AtmosGCMConfigType,
    GridInterpolated,
    "kg m^-3",
    "air density",
    "air_density",
) do (atmos::AtmosModel, states, curr_time)
    states.prognostic.ρ
end

@pointwise_diagnostic(
    "temp",
    AtmosGCMConfigType,
    GridInterpolated,
    "K",
    "air temperature",
    "air_temperature",
) do (atmos::AtmosModel, states, curr_time)
    states.thermodynamic.temp
end

@pointwise_diagnostic(
    "pres",
    AtmosGCMConfigType,
    GridInterpolated,
    "Pa",
    "air pressure",
    "air_pressure",
) do (atmos::AtmosModel, states, curr_time)
    states.thermodynamic.pres
end

@pointwise_diagnostic(
    "thd",
    AtmosGCMConfigType,
    GridInterpolated,
    "K",
    "dry potential temperature",
    "air_potential_temperature",
) do (atmos::AtmosModel, states, curr_time)
    states.thermodynamic.θ_dry
end

@pointwise_diagnostic(
    "et",
    AtmosGCMConfigType,
    GridInterpolated,
    "J kg^-1",
    "total specific energy",
    "specific_dry_energy_of_air",
) do (atmos::AtmosModel, states, curr_time)
    states.prognostic.ρe / states.prognostic.ρ
end

@pointwise_diagnostic(
    "ei",
    AtmosGCMConfigType,
    GridInterpolated,
    "J kg^-1",
    "specific internal energy",
    "internal_energy",
) do (atmos::AtmosModel, states, curr_time)
    states.thermodynamic.e_int
end

@pointwise_diagnostic(
    "ht",
    AtmosGCMConfigType,
    GridInterpolated,
    "J kg^-1",
    "specific enthalpy based on total energy",
    "",
) do (atmos::AtmosModel, states, curr_time)
    states.thermodynamic.h_tot
end

@pointwise_diagnostic(
    "hi",
    AtmosGCMConfigType,
    GridInterpolated,
    "J kg^-1",
    "specific enthalpy based on internal energy",
    "atmosphere_enthalpy_content",
) do (atmos::AtmosModel, states, curr_time)
    states.thermodynamic.h_int
end

#= TODO
@XXX_diagnostic(
    "vort",
    AtmosGCMConfigType,
    GridInterpolated,
    "s^-1",
    "vertical component of relative velocity",
    "atmosphere_relative_velocity",
) do (atmos::AtmosModel, states, curr_time)
end
=#

@pointwise_diagnostic(
    "qt",
    AtmosGCMConfigType,
    GridInterpolated,
    "kg kg^-1",
    "mass fraction of total water in air (qv+ql+qi)",
    "mass_fraction_of_water_in_air",
) do (atmos::AtmosModelm::Union{EquilMoist, NonEquilMoist}, states, curr_time)
    states.prognostic.moisture.ρq_tot / states.prognostic.ρ
end

@pointwise_diagnostic(
    "ql",
    AtmosGCMConfigType,
    GridInterpolated,
    "kg kg^-1",
    "mass fraction of liquid water in air",
    "mass_fraction_of_cloud_liquid_water_in_air",
) do (atmos::AtmosModelm::Union{EquilMoist, NonEquilMoist}, states, curr_time)
    states.thermodynamic.moisture.q_liq
end

@pointwise_diagnostic(
    "qv",
    AtmosGCMConfigType,
    GridInterpolated,
    "kg kg^-1",
    "mass fraction of water vapor in air",
    "specific_humidity",
) do (atmos::AtmosModelm::Union{EquilMoist, NonEquilMoist}, states, curr_time)
    states.thermodynamic.moisture.q_vap
end

@pointwise_diagnostic(
    "qi",
    AtmosGCMConfigType,
    GridInterpolated,
    "kg kg^-1",
    "mass fraction of ice in air",
    "mass_fraction_of_cloud_ice_in_air",
) do (atmos::AtmosModelm::Union{EquilMoist, NonEquilMoist}, states, curr_time)
    states.thermodynamic.moisture.q_ice
end

@pointwise_diagnostic(
    "thv",
    AtmosGCMConfigType,
    GridInterpolated,
    "K",
    "virtual potential temperature",
    "virtual_potential_temperature",
) do (atmos::AtmosModelm::Union{EquilMoist, NonEquilMoist}, states, curr_time)
    states.thermodynamic.moisture.θ_vir
end

@pointwise_diagnostic(
    "thl",
    AtmosGCMConfigType,
    GridInterpolated,
    "K",
    "liquid-ice potential temperature",
    "",
) do (atmos::AtmosModelm::Union{EquilMoist, NonEquilMoist}, states, curr_time)
    states.thermodynamic.moisture.θ_liq_ice
end
