#!/usr/bin/env julia --project
using ArgParse

using ClimateMachine
using ClimateMachine.Atmos
using ClimateMachine.Orientations
using ClimateMachine.ConfigTypes
using ClimateMachine.Diagnostics
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.TemperatureProfiles
using ClimateMachine.Thermodynamics
using ClimateMachine.TurbulenceClosures
using ClimateMachine.VariableTemplates
using StaticArrays
using Test
using CLIMAParameters
using CLIMAParameters.Atmos.SubgridScale: C_smag
using CLIMAParameters.Planet: R_d, cp_d, cv_d, MSLP, grav
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet();

function init_risingbubble!(problem, bl, state, aux, (x, y, z), t)
    ## Problem float-type
    FT = eltype(state)

    ## Unpack constant parameters
    R_gas::FT = R_d(bl.param_set)
    c_p::FT = cp_d(bl.param_set)
    c_v::FT = cv_d(bl.param_set)
    p0::FT = MSLP(bl.param_set)
    _grav::FT = grav(bl.param_set)
    γ::FT = c_p / c_v

    ## Define bubble center and background potential temperature
    xc::FT = 5000
    yc::FT = 1000
    zc::FT = 2000
    r = sqrt((x - xc)^2 + (z - zc)^2)
    rc::FT = 2000
    θamplitude::FT = 2
    q_tot_amplitude::FT = 1e-3

    ## TODO: clean this up, or add convenience function:
    ## This is configured in the reference hydrostatic state
    θ_ref::FT = bl.ref_state.virtual_temperature_profile.T_surface

    ## Add the thermal perturbation:
    Δθ::FT = 0
    Δq_tot::FT = 0
    if r <= rc
        Δθ = θamplitude * (1.0 - r / rc)
        Δq_tot = q_tot_amplitude * (1.0 - r / rc)
    end

    ## Compute perturbed thermodynamic state:
    θ = θ_ref + Δθ                                      # potential temperature
    q_pt = PhasePartition(Δq_tot)
    R_m = gas_constant_air(bl.param_set, q_pt)
    _cp_m = cp_m(bl.param_set, q_pt)
    _cv_m = cv_m(bl.param_set, q_pt)
    π_exner = FT(1) - _grav / (_cp_m * θ) * z           # exner pressure
    ρ = p0 / (R_gas * θ) * (π_exner)^(_cv_m / R_m)      # density
    T = θ * π_exner

    if bl.moisture isa EquilMoist
        e_int = internal_energy(bl.param_set, T, q_pt)
        ts = PhaseEquil(bl.param_set, e_int, ρ, Δq_tot)
    else
        e_int = internal_energy(bl.param_set, T)
        ts = PhaseDry(bl.param_set, e_int, ρ)
    end
    ρu = SVector(FT(0), FT(0), FT(0))                   # momentum
    ## State (prognostic) variable assignment
    e_kin = FT(0)                                       # kinetic energy
    e_pot = gravitational_potential(bl, aux)            # potential energy
    ρe_tot = ρ * total_energy(e_kin, e_pot, ts)         # total energy
    ρq_tot = ρ * Δq_tot                                 # total water specific humidity

    ## Assign State Variables
    state.ρ = ρ
    state.ρu = ρu
    state.ρe = ρe_tot
    if bl.moisture isa EquilMoist
        state.moisture.ρq_tot = ρq_tot
    end
end

function config_risingbubble(FT, N, resolution, xmax, ymax, zmax, with_moisture)

    ode_solver = ClimateMachine.ExplicitSolverType(
        solver_method = LSRK144NiegemannDiehlBusch,
    )

    T_surface = FT(300)
    T_min_ref = FT(0)
    T_profile = DryAdiabaticProfile{FT}(param_set, T_surface, T_min_ref)
    ref_state = HydrostaticState(T_profile)

    _C_smag = FT(C_smag(param_set))

    #if with_moisture
    moisture = EquilMoist{FT}()
    #else
    #    moisture = DryModel()
    #end
    model = AtmosModel{FT}(
        AtmosLESConfigType,
        param_set;
        init_state_prognostic = init_risingbubble!,
        ref_state = ref_state,
        turbulence = SmagorinskyLilly(_C_smag),
        moisture = moisture,
        source = (Gravity(),),
        tracers = NoTracers(),
    )

    config = ClimateMachine.AtmosLESConfiguration(
        "RisingBubble",
        N,
        resolution,
        xmax,
        ymax,
        zmax,
        param_set,
        init_risingbubble!,
        solver_type = ode_solver,
        model = model,
    )
    return config
end

function config_diagnostics(driver_config)
    FT = Float64
    interval = "50ssecs"
    boundaries = [
        FT(0.0) FT(0.0) FT(0.0)
        FT(10000) FT(500) FT(10000)
    ]
    resolution = (FT(100), FT(100), FT(100))
    interpol = ClimateMachine.InterpolationConfiguration(
        driver_config,
        boundaries,
        resolution,
    )
    dgngrp = setup_atmos_default_diagnostics(
        AtmosLESConfigType(),
        interval,
        driver_config.name,
    )
    state_dgngrp = setup_dump_state_diagnostics(
        AtmosLESConfigType(),
        interval,
        driver_config.name,
        interpol = interpol,
    )
    aux_dgngrp = setup_dump_aux_diagnostics(
        AtmosLESConfigType(),
        interval,
        driver_config.name,
        interpol = interpol,
    )
    return ClimateMachine.DiagnosticsConfiguration([
        dgngrp,
        state_dgngrp,
        aux_dgngrp,
    ])
end

function main()
    # add a command line argument to specify whether to use a moist setup
    rb_args = ArgParseSettings(autofix_names = true)
    add_arg_group!(rb_args, "RisingBubble")
    @add_arg_table! rb_args begin
        "--with-moisture"
        help = "use a moist setup"
        action = :store_const
        constant = true
        default = false
    end
    cl_args = ClimateMachine.init(parse_clargs = true, custom_clargs = rb_args)
    with_moisture = cl_args["with_moisture"]

    FT = Float64
    N = 4
    Δh = FT(125)
    Δv = FT(125)
    resolution = (Δh, Δh, Δv)
    xmax = FT(10000)
    ymax = FT(500)
    zmax = FT(10000)
    t0 = FT(0)
    timeend = FT(10) #FT(1000)

    CFL = FT(1.7)

    # config mystery
    driver_config =
        config_risingbubble(FT, N, resolution, xmax, ymax, zmax, with_moisture)
    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        init_on_cpu = true,
        Courant_number = CFL,
    )
    dgn_config = config_diagnostics(driver_config)
    model = driver_config.bl

    # label for ρq_tot prognostic var
    ρq_tot_idx = varsindices(
        vars_state(dg.balance_law, Prognostic(), FT),
        "moisture.ρq_tot",
    )

    # TMAR filter (without MPP)
    cb_tmar = EveryXSimulationSteps(1) do
        if odesolver isa MPPSolver
            Filters.apply!(
                state_prognostic,
                ("moisture.ρq_tot",),
                dg.grid,
                Filters.TMARFilter(),
            )
        end
        nothing
    end

    # vtk output
    output_dir = "output"
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        mkpath(output_dir)
    end
    cb_vtk = Callbacks.vtk("100steps", solver_config, output_dir, 0)

    # chatter in terminal
    cb_updates = Callbacks.show_updates("60secs", solver_config, () -> nothing)


    max_ρq_tot_init = maximum(solver_config.dg.state_auxiliary[:, ρq_tot_idx, :])
    min_ρq_tot_init = maximum(solver_config.dg.state_auxiliary[:, ρq_tot_idx, :])
    ∫ρq_tot_init = weightedsum(state_prognostic, ρq_tot_idx)

    result = ClimateMachine.invoke!(
        solver_config;
        diagnostics_config = dgn_config,
        user_callbacks = (cb_tmar_filter, cb_vtk, cb_updates()),
        check_euclidean_distance = true,
    )

    max_ρq_tot_fini = maximum(solver_config.dg.state_auxiliary[:, ρq_tot_idx, :])
    min_ρq_tot_fini = maximum(solver_config.dg.state_auxiliary[:, ρq_tot_idx, :])
    ∫ρq_tot_fini = weightedsum(state_prognostic, ρq_tot_idx)

    @info("  ")
    @info(max_ρq_tot_init, max_ρq_tot_fini)
    @info(min_ρq_tot_init, min_ρq_tot_fini)

    @info(∫ρq_tot_init, ∫ρq_tot_finit)

    #@test isapprox(result, FT(1); atol = 1.5e-3)
end

main()
