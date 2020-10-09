@intermediate_value(AtmosLESConfigType, "D_t") do (
    atmos::AtmosModel,
    states::States,
    curr_time,
    intermediates,
)
    _, D_t, _ = turbulence_tensors(
        atmos,
        states.prognostic,
        states.gradient_flux,
        states.auxiliary,
        curr_time,
    )
end

@horizontal_average(AtmosLESConfigType, "u", "m s^-1", "x-velocity", "")
@horizontal_average(AtmosLESConfigType, "v", "m s^-1", "y-velocity", "")
@horizontal_average(AtmosLESConfigType, "w", "m s^-1", "z-velocity", "")

@horizontal_average(AtmosLESConfigType, "uu")
@horizontal_average(AtmosLESConfigType, "vv")
@horizontal_average(AtmosLESConfigType, "ww")
@horizontal_average(AtmosLESConfigType, "www")
@horizontal_average(AtmosLESConfigType, "wu")
@horizontal_average(AtmosLESConfigType, "wv")

@horizontal_average_impl(
    AtmosLESConfigType,
    "u",
    "v",
    "w",
    "uu",
    "vv",
    "ww",
    "www",
    "wu",
    "wv",
) do (atmos::AtmosModel, states::States, curr_time, intermediates)
    u = states.prognostic.ρu[1]
    v = states.prognostic.ρu[2]
    w = states.prognostic.ρu[3]
    uu = u^2 / states.prognostic.ρ
    vv = v^2 / states.prognostic.ρ
    ww = w^2 / states.prognostic.ρ
    www = w^3 / states.prognostic.ρ^2
    wu = w * u / states.prognostic.ρ
    wv = w * v / states.prognostic.ρ
end

@horizontal_average(AtmosLESConfigType, "avg_rho", "kg m^-3", "air density", "air_density")
@horizontal_average(AtmosLESConfigType, "rho", "kg m^-3", "air density", "air_density")

@horizontal_average(AtmosLESConfigType, "wrho")

@horizontal_average_impl(
    AtmosLESConfigType,
    "avg_rho",
    "rho",
    "wrho",
) do (
    atmos::AtmosModel,
    states::States,
    curr_time,
    intermediates,
)
    avg_rho = states.prognostic.ρ
    rho = states.prognostic.ρ * states.prognostic.ρ
    wrho = states.prognostic.ρu[3] * states.prognostic.ρ
end

@horizontal_average(AtmosLESConfigType, "temp", "K", "air temperature", "air_temperature")
@horizontal_average(AtmosLESConfigType, "pres", "Pa", "air pressure", "air_pressure")
@horizontal_average(AtmosLESConfigType, "thd", "K", "dry potential temperature", "air_potential_temperature")
@horizontal_average(AtmosLESConfigType, "et", "J kg^-1", "total specific energy", "specific_dry_energy_of_air")
@horizontal_average(AtmosLESConfigType, "ei", "J kg^-1", "specific internal energy", "internal_energy")
@horizontal_average(AtmosLESConfigType, "ht", "J kg^-1", "specific enthalpy based on total energy", "")
@horizontal_average(AtmosLESConfigType, "hi", "J kg^-1", "specific enthalpy based on internal energy", "atmosphere_enthalpy_content")
@horizontal_average(AtmosLESConfigType, "w_ht_sgs", "kg kg^-1 m s^-1", "vertical sgs flux of total specific enthalpy", "")

@horizontal_average(AtmosLESConfigType, "eiei")
@horizontal_average(AtmosLESConfigType, "wthd")
@horizontal_average(AtmosLESConfigType, "wei")

@horizontal_average_impl(
    AtmosLESConfigType,
    "temp",
    "pres",
    "et",
    "ei",
    "ht",
    "hi",
    "w_ht_sgs",
    "eiei",
    "wthd",
    "wei",
) do (
    atmos::AtmosModel,
    states::States,
    curr_time,
    intermediates,
)
    temp = intermediates.temp * states.prognostic.ρ
    pres = intermediates.pres * states.prognostic.ρ
    thd = intermediates.θ_dry * states.prognostic.ρ
    et = states.prognostic.ρe
    ei = intermediates.e_int * states.prognostic.ρ
    ht = intermediates.h_tot * states.prognostic.ρ
    hi = intermediates.h_int * states.prognostic.ρ
    d_h_tot = -(intermediates.D_t) .* states.gradient_flux.∇h_tot
    w_ht_sgs = d_h_tot[end] * state.prognostic.ρ
    eiei = ei * intermediates.e_int
    wthd = states.prognostic.ρu[3] * intermediates.θ_dry
    wei = states.prognostic.ρu[3] * intermediates.e_int
end

@horizontal_average(AtmosLESConfigType, "qt", "kg kg^-1", "mass fraction of total water in air (qv+ql+qi)", "mass_fraction_of_water_in_air")
@horizontal_average(AtmosLESConfigType, "ql", "kg kg^-1", "mass fraction of liquid water in air", "mass_fraction_of_cloud_liquid_water_in_air")
@horizontal_average(AtmosLESConfigType, "qi", "kg kg^-1", "mass fraction of ice in air", "mass_fraction_of_cloud_ice_in_air")
@horizontal_average(AtmosLESConfigType, "qv", "kg kg^-1", "mass fraction of water vapor in air", "specific_humidity")
@horizontal_average(AtmosLESConfigType, "thv", "K", "virtual potential temperature", "virtual_potential_temperature")
@horizontal_average(AtmosLESConfigType, "thl", "K", "liquid-ice potential temperature", "")
@horizontal_average(AtmosLESConfigType, "w_qt_sgs", "kg kg^-1 m s^-1", "vertical sgs flux of total specific humidity", "")

@horizontal_average(AtmosLESConfigType, "qtqt")
@horizontal_average(AtmosLESConfigType, "thlthl")
@horizontal_average(AtmosLESConfigType, "wqt")
@horizontal_average(AtmosLESConfigType, "wql")
@horizontal_average(AtmosLESConfigType, "wqi")
@horizontal_average(AtmosLESConfigType, "wqv")
@horizontal_average(AtmosLESConfigType, "wthv")
@horizontal_average(AtmosLESConfigType, "wthl")
@horizontal_average(AtmosLESConfigType, "qtthl")
@horizontal_average(AtmosLESConfigType, "qtei")

@horizontal_average_impl(
    AtmosLESConfigType,
    "qt",
    "ql",
    "qi",
    "qv",
    "thv",
    "thl",
    "w_qt_sgs",
    "qtqt",
    "thlthl",
    "wqt",
    "wql",
    "wqi",
    "wqv",
    "wthv",
    "wthl",
    "qtthl",
    "qtei",
) do (
    m::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    intermediates
)
    qt = states.prognostic.moisture.ρq_tot
    ql = intermediates.moisture.q_liq * states.prognostic.ρ # TODO: `moisture.`?
    qi = intermediates.moisture.q_ice * states.prognostic.ρ
    qv = intermediates.moisture.q_vap * states.prognostic.ρ
    thv = intermediates.moisture.θ_vir * states.prognostic.ρ
    thl = intermediates.moisture.θ_liq_ice * states.prognostic.ρ
    d_q_tot = (-intermediates.D_t) .* states.gradient_flux.moisture.∇q_tot
    w_qt_sgs = d_q_tot[end] * states.prognostic.ρ
    qtqt = qt * (states.prognostic.moisture.ρq_tot / states.prognostic.ρ)
    thlthl = thl * intermediates.moisture.θ_liq_ice
    w = states.prognostic.ρu[3] / states.prognostic.ρ
    wqt = states.prognostic.ρu[3] * qt / states.prognostic.ρ
    wql = states.prognostic.ρu[3] * intermediates.moisture.q_liq
    wqi = states.prognostic.ρu[3] * intermediates.moisture.q_ice
    wqv = states.prognostic.ρu[3] * intermediates.moisture.q_vap
    wthv = states.prognostic.ρu[3] * intermediates.moisture.θ_vir
    wthl = states.prognostic.ρu[3] * intermediates.moisture.θ_liq_ice
    qtthl = qt * intermediates.moisture.θ_liq_ice
    qtei = qt * intermediates.e_int
end

#= TODO
    Variables["cld_frac"] = DiagnosticVariable(
        "cld_frac",
        diagnostic_var_attrib(
            "",
            "cloud fraction",
            "cloud_area_fraction_in_atmosphere_layer",
        ),
    )
    Variables["cld_cover"] = DiagnosticVariable(
        "cld_cover",
        diagnostic_var_attrib("", "cloud cover", "cloud_area_fraction"),
    )
    Variables["cld_top"] = DiagnosticVariable(
        "cld_top",
        diagnostic_var_attrib("m", "cloud top", "cloud_top_altitude"),
    )
    Variables["cld_base"] = DiagnosticVariable(
        "cld_base",
        diagnostic_var_attrib("m", "cloud base", "cloud_base_altitude"),
    )
    Variables["lwp"] = DiagnosticVariable(
        "lwp",
        diagnostic_var_attrib(
            "kg m^-2",
            "liquid water path",
            "atmosphere_mass_content_of_cloud_condensed_water",
        ),
    )
    Variables["core_frac"] = DiagnosticVariable(
        "core_frac",
        diagnostic_var_attrib("", "cloud core fraction", ""),
    )
    Variables["u_core"] = DiagnosticVariable(
        "u_core",
        diagnostic_var_attrib("m s^-1", "cloud core x-velocity", ""),
    )
    Variables["v_core"] = DiagnosticVariable(
        "v_core",
        diagnostic_var_attrib("m s^-1", "cloud core y-velocity", ""),
    )
    Variables["w_core"] = DiagnosticVariable(
        "w_core",
        diagnostic_var_attrib("m s^-1", "cloud core z-velocity", ""),
    )
    Variables["avg_rho_core"] = DiagnosticVariable(
        "avg_rho_core",
        diagnostic_var_attrib("kg m^-3", "cloud core air density", ""),
    )
    Variables["rho_core"] = DiagnosticVariable(
        "rho_core",
        diagnostic_var_attrib("kg m^-3", "cloud core (density-averaged) air density", ""),
    )
    Variables["qt_core"] = DiagnosticVariable(
        "qt_core",
        diagnostic_var_attrib("kg m^-3", "cloud core total specific humidity", ""),
    )
    Variables["ql_core"] = DiagnosticVariable(
        "ql_core",
        diagnostic_var_attrib("kg m^-3", "cloud core liquid water specific humidity", ""),
    )
    Variables["thv_core"] = DiagnosticVariable(
        "thv_core",
        diagnostic_var_attrib("K", "cloud core virtual potential temperature", ""),
    )
    Variables["thl_core"] = DiagnosticVariable(
        "thl_core",
        diagnostic_var_attrib("K", "cloud core liquid-ice potential temperature", ""),
    )
    Variables["ei_core"] = DiagnosticVariable(
        "ei_core",
        diagnostic_var_attrib("J kg-1", "cloud core specific internal energy", ""),
    )
    Variables["var_u_core"] = DiagnosticVariable(
        "var_u_core",
        diagnostic_var_attrib("m^2 s^-2", "cloud core variance of x-velocity", ""),
    )
    Variables["var_v_core"] = DiagnosticVariable(
        "var_v_core",
        diagnostic_var_attrib("m^2 s^-2", "cloud core variance of y-velocity", ""),
    )
    Variables["var_w_core"] = DiagnosticVariable(
        "var_w_core",
        diagnostic_var_attrib("m^2 s^-2", "cloud core variance of z-velocity", ""),
    )
    Variables["var_qt_core"] = DiagnosticVariable(
        "var_qt_core",
        diagnostic_var_attrib(
            "kg^2 kg^-2",
            "cloud core variance of total specific humidity",
            "",
        ),
    )
    Variables["var_thl_core"] = DiagnosticVariable(
        "var_thl_core",
        diagnostic_var_attrib(
            "K^2",
            "cloud core variance of liquid-ice potential temperature",
            "",
        ),
    )
    Variables["var_ei_core"] = DiagnosticVariable(
        "var_ei_core",
        diagnostic_var_attrib(
            "J^2 kg^-2",
            "cloud core variance of specific internal energy",
            "",
        ),
    )
    Variables["cov_w_rho_core"] = DiagnosticVariable(
        "cov_w_rho_core",
        diagnostic_var_attrib(
            "kg m^-2 s^-1",
            "cloud core vertical eddy flux of density",
            "",
        ),
    )
    Variables["cov_w_qt_core"] = DiagnosticVariable(
        "cov_w_qt_core",
        diagnostic_var_attrib(
            "kg kg^-1 m s^-1",
            "cloud core vertical eddy flux of specific humidity",
            "",
        ),
    )
    Variables["cov_w_thl_core"] = DiagnosticVariable(
        "cov_w_thl_core",
        diagnostic_var_attrib(
            "K m s^-1",
            "cloud core vertical eddy flux of liquid-ice potential temperature",
            "",
        ),
    )
    Variables["cov_w_ei_core"] = DiagnosticVariable(
        "cov_w_ei_core",
        diagnostic_var_attrib(
            "J kg^-1 m^-1 s^-1",
            "cloud core vertical eddy flux of specific internal energy",
            "",
        ),
    )
    Variables["cov_qt_thl_core"] = DiagnosticVariable(
        "cov_qt_thl_core",
        diagnostic_var_attrib(
            "kg kg^-1 K",
            "cloud core covariance of total specific humidity and liquid-ice potential temperature",
            "",
        ),
    )
    Variables["cov_qt_ei_core"] = DiagnosticVariable(
        "cov_qt_ei_core",
        diagnostic_var_attrib(
            "kg kg^-1 J kg^-1",
            "cloud core covariance of total specific humidity and specific internal energy",
            "",
        ),
    )
    Variables["E_k"] = DiagnosticVariable(
        "E_k",
        diagnostic_var_attrib(
            "",
            "volumetrically-averaged dimensionless kinetic energy",
            "",
        ),
    )
    Variables["dE"] = DiagnosticVariable(
        "dE",
        diagnostic_var_attrib(
            "",
            "volumetrically-averaged kinetic energy dissipation",
            "",
        ),
    )
    Variables["mass_loss"] =
        DiagnosticVariable("mass_loss", diagnostic_var_attrib("", "", ""))
    Variables["energy_loss"] =
        DiagnosticVariable("energy_loss", diagnostic_var_attrib("", "", ""))
=#
