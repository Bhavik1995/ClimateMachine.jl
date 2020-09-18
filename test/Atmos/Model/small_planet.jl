using ClimateMachine
ClimateMachine.init()
using ClimateMachine.Atmos
using ClimateMachine.ConfigTypes
using ClimateMachine.Diagnostics
using ClimateMachine.Mesh.Interpolation
using CLIMAParameters
using CLIMAParameters.Planet: planet_radius
import CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
CLIMAParameters.Planet.planet_radius(::EarthParameterSet) = 10^4
function main()
    FT,n_horz,n_vert,poly_order = Float64,12,6,3
    model = AtmosModel{FT}(AtmosGCMConfigType, param_set; init_state_prognostic = x->x)
    domain_height = FT(10^4*4)
    _planet_radius = FT(planet_radius(param_set))
    driver_config = ClimateMachine.AtmosGCMConfiguration(
        "SmallPlanet", poly_order, (n_horz, n_vert),
        domain_height, param_set, x->x; model = model,
    )
    info = driver_config.config_info
    boundaries = [ FT(-90) FT(-180) _planet_radius
        FT(90) FT(180) FT(_planet_radius + info.domain_height)]
    resolution = (FT(2), FT(2), FT(1000)) # in (deg, deg, m)
    interpol = ClimateMachine.InterpolationConfiguration(
        driver_config,
        boundaries,
        resolution,
    )
    return nothing
end
main()
