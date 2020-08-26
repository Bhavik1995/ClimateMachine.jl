data_name(::Prognostic) = "prog"
data_name(::Auxiliary) = "aux"
data_name(::GradientFlux) = "grad_flux"

get_data(solver_config, ::Prognostic) = solver_config.Q.data
get_data(solver_config, ::Auxiliary) = solver_config.dg.state_auxiliary.data
get_data(solver_config, ::GradientFlux) = solver_config.dg.state_gradient_flux.data

using HDF5
function export_state(
    solver_config,
    output_dir;
    state_types = (Prognostic(), Auxiliary()),
)
    FT = eltype(solver_config.Q)
    mkpath(output_dir)
    for st in state_types
        data = get_data(solver_config, st)
        name = data_name(st)
        h5write(joinpath(output_dir, name*".h5")   , "mydata/$name", data)
    end
end

function import_state(
    import_dir;
    state_types = (Prognostic(), Auxiliary()),
)
    data = Dict()
    for st in state_types
        name = data_name(st)
        data[st] = h5read(joinpath(import_dir, name*".h5"), "mydata/$name")
    end
    return data
end
