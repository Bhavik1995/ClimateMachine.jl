"""
    ConfigTypes

Module containing ClimateMachine configuration types.
"""
module ConfigTypes

export ClimateMachineConfigType,
    AtmosLESConfigType,
    AtmosGCMConfigType,
    OceanBoxGCMConfigType,
    OceanSplitExplicitConfigType,
    SingleStackConfigType

abstract type ClimateMachineConfigType end
abstract type AtmosConfigType <: ClimateMachineConfigType end
struct AtmosLESConfigType <: AtmosConfigType end
struct AtmosGCMConfigType <: AtmosConfigType end
struct OceanBoxGCMConfigType <: ClimateMachineConfigType end
struct OceanSplitExplicitConfigType <: ClimateMachineConfigType end
struct SingleStackConfigType <: ClimateMachineConfigType end

end
