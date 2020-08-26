using ClimateMachine
const clima_dir = dirname(dirname(pathof(ClimateMachine)));
include(joinpath(clima_dir, "experiments", "AtmosLES", "bomex_model.jl"))

function main(::Type{FT}) where {FT}

    # DG polynomial order
    N = 1
    nelem_vert = 50

    # Prescribe domain parameters
    zmax = FT(3000)

    t0 = FT(0)

    # For a full-run, please set the timeend to 3600*6 seconds
    # For the test we set this to == 30 minutes
    # timeend = FT(13.805585)
    # timeend = FT(400)
    timeend = FT(400)
    #timeend = FT(3600 * 6)
    CFLmax = FT(0.90)
    FT = Float64
    config_type = SingleStackConfigType

    # Choose default IMEX solver
    ode_solver_type = ClimateMachine.IMEXSolverType()

    N_updrafts = 1
    N_quad = 3
    turbconv = EDMF(FT, N_updrafts, N_quad)

    model = bomex_model(FT, config_type, zmax, surface_flux, turbconv)

    # Assemble configuration
    driver_config = ClimateMachine.SingleStackConfiguration(
        "BOMEX_SINGLE_STACK",
        N,
        nelem_vert,
        zmax,
        param_set,
        model;
        solver_type = ode_solver_type,
    )

    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        init_on_cpu = true,
        Courant_number = CFLmax,
    )

    dgn_config = config_diagnostics(driver_config)

    cbtmarfilter = GenericCallbacks.EveryXSimulationSteps(1) do
        Filters.apply!(
            solver_config.Q,
            (
                "moisture.ρq_tot",
                "turbconv.environment.ρatke",
                "turbconv.environment.ρaθ_liq_cv",
                "turbconv.environment.ρaq_tot_cv",
                "turbconv.updraft",
            ),
            solver_config.dg.grid,
            TMARFilter(),
        )
        nothing
    end

    # State variable
    Q = solver_config.Q
    # Volume geometry information
    vgeo = driver_config.grid.vgeo
    M = vgeo[:, Grids._M, :]
    # Unpack prognostic vars
    ρ₀ = Q.ρ
    ρe₀ = Q.ρe
    # DG variable sums
    Σρ₀ = sum(ρ₀ .* M)
    Σρe₀ = sum(ρe₀ .* M)

    # -------------------------- Quick & dirty diagnostics. TODO: replace with proper diagnostics

    grid = driver_config.grid
    output_dir = ClimateMachine.Settings.output_dir
    @show output_dir
    state_types = (Prognostic(), Auxiliary(), GradientFlux())
    all_data = [dict_of_nodal_states(solver_config, ["z"], state_types)]
    time_data = FT[0]

    export_state_plots(
        solver_config,
        all_data,
        time_data,
        joinpath(clima_dir, "output", "ICs");
        state_types = state_types,
    )

    # Define the number of outputs from `t0` to `timeend`
    n_outputs = 4
    # This equates to exports every ceil(Int, timeend/n_outputs) time-step:
    every_x_simulation_time = ceil(Int, timeend / n_outputs)

    cb_data_vs_time =
        GenericCallbacks.EveryXSimulationTime(every_x_simulation_time) do
            # cb_data_vs_time = GenericCallbacks.EveryXSimulationSteps(1) do
            push!(
                all_data,
                dict_of_nodal_states(solver_config, ["z"], state_types),
            )
            push!(time_data, gettime(solver_config.solver))
            @show gettime(solver_config.solver)
            @show getsteps(solver_config.solver)
            nothing
        end
    # --------------------------

    cb_check_cons = GenericCallbacks.EveryXSimulationSteps(3000) do
        Q = solver_config.Q
        δρ = (sum(Q.ρ .* M) - Σρ₀) / Σρ₀
        δρe = (sum(Q.ρe .* M) .- Σρe₀) ./ Σρe₀
        @show (abs(δρ))
        @show (abs(δρe))
        @test (abs(δρ) <= 0.001)
        @test (abs(δρe) <= 0.0025)
        nothing
    end

    result = ClimateMachine.invoke!(
        solver_config;
        diagnostics_config = dgn_config,
        user_callbacks = (cbtmarfilter, cb_check_cons, cb_data_vs_time),
        check_euclidean_distance = true,
        # numberofsteps = 291
    )
    push!(all_data, dict_of_nodal_states(solver_config, ["z"], state_types))
    push!(time_data, gettime(solver_config.solver))

    export_state_plots(
        solver_config,
        all_data,
        time_data,
        joinpath(clima_dir, "output", "runtime");
        state_types = state_types,
    )

    @test !isnan(norm(Q))
    return solver_config, all_data, time_data
end

solver_config, all_data, time_data = main(Float64)

accept_new_solution = false
include("io_state.jl")
if accept_new_solution
    export_state(solver_config, "bomex_edmf_output")
end

bomex_data = import_state("bomex_edmf_output")
using Test
@test all(solver_config.Q.data .≈ bomex_data[Prognostic()])
# @test all(solver_config.dg.state_auxiliary.data .≈ bomex_data[Auxiliary()]) # Has some NaNs

