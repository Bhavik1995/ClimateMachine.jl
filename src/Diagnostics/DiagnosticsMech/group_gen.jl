# Generate `setup_$(name)(...)` to create the `DiagnosticsGroup`
# for $name.
function generate_setup(name, config_type, params_type, dvars)
    setup_name = Symbol("setup_", name)
    init_name = Symbol(name, "_init")
    collect_name = Symbol(name, "_collect")
    fini_name = Symbol(name, "_fini")
    quote
        function $(setup_name)(
            ::$config_type,
            params::$params_type,
            interval::String,
            out_prefix::String,
            writer = NetCDFWriter(),
            interpol = nothing,
        ) where {
            $config_type <: ClimateMachineConfigType,
            $params_type <: Union{Nothing, DiagnosticsGroupParams},
        }
            return DiagnosticsGroup(
                $(name),
                Diagnostics.$(init_name),
                Diagnostics.$(collect_name),
                Diagnostics.$(fini_name),
                interval,
                out_prefix,
                writer,
                interpol,
                $(any(v -> v.layout isa LayoutDG, dvars)),
                params,
            )
        end
    end
end

# Generate the `dims` dictionary for `Writers.init_data`.
function generate_init_dims(name, config_type, dvars)
    # Set up an error for when no InterpolationTopology is specified but there's
    # a diagnostic variable with a Cartesian layout.
    err_ex = quote end
    if any(v -> v.layout ∉ (LayoutScalar, LayoutDG), dvars)
        err_ex = quote
            throw(
                ArgumentError(
                    "Diagnostic variables with LayoutCartesian " *
                    "require an InterpolationTopology"
                )
            )
        end
    end

    # Add a `z` dimension for horizontal averages.
    add_z_dim_ex = quote end
    if any(v -> v.kind isa HorizontalAverage, dvars)
        add_z_dim_ex = quote
            dims["z"] = (AtmosCollected.zvals, Dict())
        end
    end

    quote
        dims = dimensions(interpol)
        if isempty(dims)
            $(err_ex)
            $(add_z_dim_ex)
        elseif interpol isa InterpolationCubedSphere
            # Adjust `level` on the sphere.
            level_val = dims["level"]
            dims["level"] = (
                level_val[1] .- FT(planet_radius(Settings.param_set)),
                level_val[2],
            )
        end
        dims
    end
end

# Generate the `vars` dictionary for `Writers.init_data`.
function generate_init_vars(name, config_type, dvars)
    elems = ()

    for dvar in dvars
        elems = (
            elems...,
            dv_name($config_type, $dvar) => (
                dv_dims($config_type, $dvar),
                FT,
                dv_attrib($config_type, $dvar),
            )
        )
    end

    quote
        OrderedDict($(elems...))
    end
end

# Generate `Diagnostics.$(name)_init(...)` to initialize the
# `DiagnosticsGroup`.
function generate_init(name, config_type, dvars)
    init_name = Symbol(name, "_init")
    quote
        function $(esc(init_name))(dgngrp, curr_time)
            mpicomm = Settings.mpicomm
            mpirank = MPI.Comm_rank(mpicomm)
            dg = Settings.dg
            bl = dg.balance_law
            Q = Settings.Q
            FT = eltype(Q)
            interpol = dgngrp.interpol

            if dgngrp.onetime
                atmos_collect_onetime(Settings.mpicomm, Settings.dg, Settings.Q)
            end

            if mpirank == 0
                dims = $(generate_init_dims(name, config_type, dvars))
                vars = $(generate_init_vars(name, config_type, dvars))

                # create the output file
                dprefix = @sprintf(
                    "%s_%s_%s",
                    dgngrp.out_prefix,
                    dgngrp.name,
                    Settings.starttime,
                )
                dfilename = joinpath(Settings.output_dir, dprefix)
                init_data(dgngrp.writer, dfilename, dims, vars)
            end

            return nothing
        end
    end
end

# Generate `Diagnostics.$(name)_collect(...)`, which performs a
# collection of all the diagnostic variables in the group and writes
# them out.
function generate_collect(name, vars)
    collect_name = Symbol(name, "_collect")
    quote
        function $(esc(collect_name))(dgngrp, curr_time)
            if mpirank == 0
                varvals = OrderedDict()
                # ...
                append_data(dgngrp.writer, varvals, curr_time)
            end

            return nothing
        end
    end
end

# Generate `Diagnostics.$(name)_fini(...)`, which does nothing
# right now.
function generate_fini(name, vars)
    fini_name = Symbol(name, "_fini")
    quote
        function $(esc(fini_name))(dgngrp, curr_time) end
    end
end
