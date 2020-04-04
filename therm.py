#!/usr/bin/env python3

"""
Flow / thermal calculations to design thermal desorption olfactometer.
"""

import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI
import pint


ureg = pint.UnitRegistry()

safety_factor = 1.2

ambient_temp_c = 25.0
ambient_pressure = (1 * ureg.atm).to('pascal')

# Highest recommended temperature (some sorbents had lower values) from some
# (Restek? Gerstel?) sorpent selection guide PDF:
# (I think it belonged to one of the Carbotrap variants)
desorption_temp_c = 350.0

# From Sigma-Aldrich absorbent selection guide
# TODO unclear what ID was in these plots, but presumably in 1-5mm range
# in table on page 23 (may be able to guess from values below @ 3,4mm)
#adsorbent_mesh_size2resistance = {
#    '80-100': 

# "Carbotrap" seems to specifically refer to the 20-40 mesh size particle
# variety, from Sigma's website.
# https://www.sisweb.com/index/referenc/carboprs.htm
# Values are "back" pressures in units of "inches of water"
carbotrap_3mm = pd.DataFrame(
    index=[50, 100, 150, 200],
    columns=[25, 50, 75, 100, 125, 150, 175, 200], data=[
    [1.5, 3,   4,   6., 7.5,  9,    10.5, 12],
    [1.5, 3.5, 4.5, 7,  8,    10,   11.5, 13.5],
    [2.5, 5,   7.5, 10, 12.5, 15.5, 18,   21],
    [7.5, 15,  26,  36, None, None, None, None]
])
carbotrap_3mm.columns.name = 'ml/min'
carbotrap_3mm.index.name = 'mg'

carbotrap_4mm = pd.DataFrame(
    index=[50, 100, 150, 200, 250, 300, 350],
    columns=[25, 50, 75, 100, 125, 150, 175, 200], data=[
    [1,   2,   3,   4,   5,    5.5,  6.5,  7.5],
    [1,   2,   3,   4,   5,    5.5,  6.5,  7.5],
    [1,   2,   3,   4,   5,    5.5,  6.5,  8],
    [1,   2.5, 3.5, 5,   6.5,  7,    8,    9.5],
    [1,   2.5, 3.5, 5,   6.5,  7,    8,    9.5],
    [1.2, 2.5, 3.5, 5,   6.5,  7,    8,    10],
    [2,   4,   6,   8.5, 10.5, 13,   15,   17.5]
])
carbotrap_4mm.columns.name = 'ml/min'
carbotrap_4mm.index.name = 'mg'

# TODO why did i pick this value again? might want to move it higher, to maybe
# get more stuff delivered to the fly (just b/c hopefully better trapping
# efficiency)?
carbotrap_mass = 350 * ureg.milligram


def is_dimensionless(q):
    return q.to('').magnitude == q.magnitude


def cooling():
    # carbotrap_mass defined above so it can be shared with calculations
    # on energy requirements to heat trap
    volumetric_flowrate = 200 * ureg.milliliter / ureg.minute
    trap_id = (4.0 * ureg.millimeter).to('m')
    cooling_tube_id = (4.0 * ureg.millimeter).to('m')
    cooling_tube_temp = ureg.Quantity(ambient_temp_c, ureg.celsius).to('kelvin')

    # From "Inch of water" Wikipedia page (water @ 60F)
    inches_of_water2pascals = 248.84

    # Ambient pressure (at exit of olfactometer)
    p2_pa = ambient_pressure

    # TODO approx / measure pressure required upstream of carbotrap to produce
    # flows we want, and use that for initial density calculation
    p1 = carbotrap_4mm.loc[
        carbotrap_mass.magnitude,
        volumetric_flowrate.magnitude
    ]
    # TODO i guess i would need to add this "back" pressure to the ambient
    # pressure to get the pressure upstream of the trap?
    # TODO is the differential (before adding p2_pa) really so low?
    # TODO also factor in drop along whatever lengths of tubing?
    p1_pa = p1 * inches_of_water2pascals * ureg.pascal + p2_pa

    # Maximum temperature of the helium.
    # (Slightly above?) temperature for desorption.
    t1_k = ureg.Quantity(300, ureg.celsius).to('kelvin')
    # TODO calculate this temperature of the gas exiting the olfactometer,
    # given other values. this is probably the main quantity of interest.
    # TODO TODO or just set this to ambient, and solve for length...?
    t2_k = ureg.Quantity(ambient_temp_c, ureg.celsius).to('kelvin')

    gas = 'Helium'
    d1 = PropsSI('D', 'T', t1_k.magnitude, 'P', p1_pa.magnitude, gas)
    d2 = PropsSI('D', 'T', t2_k.magnitude, 'P', p2_pa.magnitude, gas)
    d1 = d1 * ureg['kg/m**3']
    d2 = d2 * ureg['kg/m**3']

    # "Dynamic" viscosity.
    v1 = PropsSI('V', 'T', t1_k.magnitude, 'P', p1_pa.magnitude, gas)
    v2 = PropsSI('V', 'T', t2_k.magnitude, 'P', p2_pa.magnitude, gas)
    v1 = v1 * ureg['Pa s']
    v2 = v2 * ureg['Pa s']

    tube_crosssection_area = np.pi * (cooling_tube_id / 2)**2
    linear_flowrate = (volumetric_flowrate / tube_crosssection_area).to('m/s')

    # Can also be defind w/ "kinematic viscosity", which is just:
    # dynamic viscosity / density.
    reynolds =  d1 * linear_flowrate * cooling_tube_id / v1
    assert is_dimensionless(reynolds)
    reynolds = reynolds.to('')

    # Thresholds taken from Wikipedia page on Reynold's number.
    # Applies to "fully developed" flow, which is flow some threshold number
    # of tube diameters into the tube.
    laminar = False
    print('Flow in cooling tube will be ', end='')
    if reynolds.magnitude < 2300:
        print('laminar', end='')
        laminar = True
    elif 2300 <= reynolds.magnitude <= 2900:
        print('partially turbulent', end='')
    else:
        print('turbulent', end='')
    print(f' (Re={reynolds.magnitude:.1f})')

    # Thermal conductivity of the fluid.
    k1 = PropsSI('CONDUCTIVITY', 'T', t1_k.magnitude, 'P', p1_pa.magnitude, gas)
    k1 = k1 * ureg['W/m/K']
    
    # TODO do i want this or CVMASS / CPOMASS?
    # (constant pressure, constant volume, ideal gas const pressure...)
    # (they all have the same units)
    cp1 = PropsSI('C', 'T', t1_k.magnitude, 'P', p1_pa.magnitude, gas)
    # first and last are basically the same
    '''
    print(cp1)
    print(PropsSI('CVMASS', 'T', t1_k.magnitude, 'P', p1_pa.magnitude, gas))
    print(PropsSI('CP0MASS', 'T', t1_k.magnitude, 'P', p1_pa.magnitude, gas))
    '''
    cp1 = cp1 * ureg['J/kg/K']

    # Formula from Wikipedia page on Prandtl number, aka Pr.
    # "dimensionless parameter representing the ratio of diffusion of
    #  momentum to diffusion of heat"
    prandtl = cp1 * v1 / k1
    assert is_dimensionless(prandtl)
    prandtl = prandtl.to('')

    # aka Re * Pr
    peclet = reynolds * prandtl
    assert is_dimensionless(peclet)
    peclet = peclet.to('')

    ub = PropsSI('V', 'T', t1_k.magnitude, 'P', p1_pa.magnitude, gas)
    uw = PropsSI('V', 'T', cooling_tube_temp.magnitude, 'P', p1_pa.magnitude,
        gas
    )

    # TODO solve for this!
    cooling_tube_length = 0.5 * ureg['m']

    # Using Sieder and Tate approximation in equation 20-27 from 5th edition
    # "Fundamentals of Momentum, Heat, and Mass Transfer" by Welty, Wicks,
    # Wilson, and Rorrer.
    # "All properties other than uw are evaluated at the bulk fluid temperature"
    # ub = viscosity of the "bulk" fluid (average temp of stream)
    # uw = viscosity of the fluid at the tube wall
    nusselt = (1.88 * (peclet * cooling_tube_id / cooling_tube_length)**(1/3) *
        (ub/uw)**0.14
    )
    assert is_dimensionless(nusselt)
    nusselt = nusselt.to('')


    # TODO TODO TODO how to factor in thermal conductivity of wall to cooling
    # calculation?
    
    # PTFE properties from Dupont Teflon PTFE properties handbook.
    # 4.6mm is listed next to this specification, although I thought it was not
    # a property that depended on length...
    # Listed in units of W/m/K
    ptfe_k = 0.25

    # Specific heat. kJ/kg/K.
    # From 20C to 260C, reported values only vary in the range [1.2, 1.5]
    ptfe_cp = 1.4

    # TODO check validity conditions for sieder-tate approx are met.
    # this source lists some, though check where these came from.
    # https://www.nuclear-power.net/nuclear-engineering/heat-transfer/convection-convective-heat-transfer/sieder-tate-equation/
    # see also http://facstaff.cbu.edu/rprice/lectures/htcoeff.html
    # ...although both of these sites list a *high* reynold's number as a
    # requirement, but my understanding from the book I used and the 1936 paper
    # is that this correlation is applicable in *laminar* flow (i.e. *low* Re)
    # it looks like there are two different sieder-tate correlations. does 
    # the laminar one even apply for tubes of the length i'm interested in?
    # one source ("Heat transfer in flow through conduits" PDF) said it only
    # applied for short tubes...

    h = nusselt * k1 / cooling_tube_id


    import ipdb; ipdb.set_trace()


# TODO and maybe convert to mm if already a pint quantity
# (in potentially other units)?
def mm(magnitude):
    return magnitude * ureg['mm']


def tube_mass(outer_diam_mm, inner_diam_mm, length_mm, density):
    """
    Length arguments are raw ints/floats in mm. `density` is a `pint`
    quantity (with associated units).
    """
    assert outer_diam_mm > inner_diam_mm
    assert density.check('[mass] [length]^-3')

    outer_diam = mm(outer_diam_mm)
    inner_diam = mm(inner_diam_mm)
    length = mm(length_mm)

    volume = np.pi * length * ((outer_diam / 2)**2  - (inner_diam / 2)**2)
    return (volume * density).to('g')


def engtoolbox_specific_heat_units(magnitude):
    # https://www.engineeringtoolbox.com/specific-heat-capacity-d_391.html
    return magnitude * ureg['J / (kg * K)']


# TODO convert to g/ml for sanity checking
def engtoolbox_density_units(magnitude):
    # https://www.engineeringtoolbox.com/density-solids-d_1265.html
    return magnitude * 10**3 * ureg['kg / m^3']


# TODO maybe it's sufficient to just heat around the sorbent trap directly,
# and have the gas only start to be heated when it enters? in which case:

# TODO TODO maybe some combination of this + gas heater calculation?
# (where maybe one of the heaters only aims for an intermediate temperature,
# and maybe the trap heater can take an arbitrarily long amount of time to reach
# that temperature, assuming it's the one at a lower temp)

# Assuming a sampling tube such as the 6mm OD x 4mm ID x 177.8mm L empty
# Gerstel tubes available from Sigma.
# (all in mm)
trap_tube_id_mm = 4
trap_tube_od_mm = 6
trap_tube_length_mm = 177.8
def min_trap_heater_power(min_power_w=None, time_s=None):
    """
    Prints a lower bound on either the heating time (`time_s`) or the required
    power (`min_power_w`), to achieve `desorption_temp_c` in the given amount of
    time.
    """
    assert min_power_w is None or time_s is None
    if min_power_w is None and time_s is None:
        time_s = 30.0
    
    if time_s is None:
        raise NotImplementedError

    steel_specific_heat = engtoolbox_specific_heat_units(490)
    steel_density = engtoolbox_density_units(7.82)

    trap_tube_material = 'glass'
    if trap_tube_material == 'glass':
        # (assuming "quartz glass") specific heat=700 J/(kg * C)
        # (sh=670 - 840 for other types of glass)
        # (from "quartz") density=2.65 10^3 kg / m^3
        trap_tube_specific_heat = engtoolbox_specific_heat_units(700)
        trap_tube_density = engtoolbox_density_units(2.65)

    elif trap_tube_material == 'steel':
        trap_tube_specific_heat = steel_specific_heat
        trap_tube_density = steel_density
    else:
        raise ValueError('invalid trap_tube_material')

    # TODO try calculation w/ both stainless and aluminum as material here
    # (some other reason to prefer one over the other here?)
    # For a metal tube that contacts the heating element, and surrounds the trap
    # Assuming any tolerance required for fit is negligible.
    trap_heater_id_mm = trap_tube_od_mm
    trap_heater_wall_thickness_mm = 3
    # Assuming this might help add fittings that seal against the trap tube.
    trap_heater_extra_length_mm = 10

    trap_heater_od_mm = trap_heater_id_mm + trap_heater_wall_thickness_mm
    trap_heater_length_mm = trap_tube_length_mm + trap_heater_extra_length_mm

    # Steel: sh=490 J/(kg * C), density=7.82 10^3 kg / m^3
    # Aluminum: sh=897 J/(kg * C), density=2.7 10^3 kg / m^3
    # (going w/ steel now b/c it may be more inert on any surfaces that 
    # might touch the sample if there is a leak + probably more importantly,
    # there may be less concern about different heat expansion if other fittings
    # are also steel [though there will probably be temperature gradients...]?)
    trap_heater_material = 'aluminum'
    if trap_heater_material == 'steel':
        trap_heater_specific_heat = steel_specific_heat
        trap_heater_density = steel_density

    elif trap_heater_material == 'aluminum':
        trap_heater_specific_heat = engtoolbox_specific_heat_units(897)
        trap_heater_density = engtoolbox_density_units(2.7)
    else:
        raise ValueError('invalid trap_heater_material')

    trap_tube_mass = tube_mass(trap_tube_od_mm, trap_tube_id_mm,
        trap_tube_length_mm, trap_tube_density
    )
    trap_heater_mass = tube_mass(trap_heater_od_mm, trap_heater_id_mm,
        trap_heater_length_mm, trap_heater_density
    )

    # carbon black (seems to be what carbotrap is) specific heat:
    # 0.165 cal / g @ 25 C
    #   https://www.reade.com/products/carbon-black
    # (assuming "graphite (carbon)") 717 J/(kg * C)
    # (pretty close to 690.36 of reade.com value if converted to same units)
    sorbent_specific_heat = engtoolbox_specific_heat_units(717)
    sorbent_mass = carbotrap_mass

    # TODO adrian recommended 50/50 this and carbotrap, right?
    # (or should they maybe be sequential, w/ stronger second, and desorption
    # flow in opposite direction...?)
    # (Tenax) (can't find specific heat for this, so just assuming the trap is
    # full of carbotrap, or that the specific heats for the two don't differ
    # much)

    # other name for this kind of quantity besides "energy_per_c"?
    trap_heater_energy_per_c = trap_heater_mass * trap_heater_specific_heat
    trap_tube_energy_per_c = trap_tube_mass * trap_tube_specific_heat
    sorbent_energy_per_c = sorbent_mass * sorbent_specific_heat

    total_energy_per_c = \
        trap_heater_energy_per_c + trap_tube_energy_per_c + sorbent_energy_per_c

    assert desorption_temp_c > ambient_temp_c
    temperature_delta = (
        ureg.Quantity(desorption_temp_c, ureg.celsius) -
        ureg.Quantity(ambient_temp_c, ureg.celsius)
    ).to('kelvin')

    # Assuming perfect efficiency, apart from safety factor.
    min_heating_energy = total_energy_per_c * temperature_delta * safety_factor

    def watts(quantity):
        """Returns magnitude of power quantity in watts"""
        return quantity.to('watt').magnitude

    time = time_s * ureg['s']
    min_heating_power_w = watts(min_heating_energy / time)
    print(f'Minimum power required to heat trap and container to '
        f'{desorption_temp_c:.0f}C in {time_s:.0f}s: '
        f'{min_heating_power_w:.0f}W'
    )
    trap_heater_power_w = watts(trap_heater_energy_per_c * temperature_delta *
        safety_factor / time
    )
    trap_tube_power_w = watts(trap_tube_energy_per_c * temperature_delta *
        safety_factor / time
    )
    sorbent_power_w = watts(sorbent_energy_per_c * temperature_delta *
        safety_factor / time
    )
    print(f' - heater body: {trap_heater_power_w:.1f}W')
    print(f' - trap tube: {trap_tube_power_w:.1f}W')
    print(f' - sorbent: {sorbent_power_w:.1f}W')

    # From High-Temperature section in:
    # https://www.mcmaster.com/heaters/heaters-for-pipes-and-tubes-5
    # These heaters are also all 
    mcmaster_wattage2length_in = {
        25: 6,
        35: 8,
        50: 12,
        100: 24,
        125: 36,
        250: 60,
        400: 96,
        500: 120
    }
    # All of the above heaters list their thickness as 0.17"
    heater_thickness = 0.17 * ureg['in']

    mcm_heater_length_in = None
    for mcm_heater_wattage in sorted(mcmaster_wattage2length_in.keys()):
        if mcm_heater_wattage >= min_heating_power_w:
            print(f'{mcm_heater_wattage}W McMaster heater is sufficient')
            mcm_heater_length_in = \
                mcmaster_wattage2length_in[mcm_heater_wattage]
            break

    if mcm_heater_length_in is None:
        print('Required power higher than any McMaster High-Temperature tube '
            'heater can provide!'
        )
        return

    # TODO TODO TODO can i PWM these heaters, or do i NEED the "variable
    # voltage" controllers? why do they sell them with plugs, otherwise?
    
    heater_body_circumference = np.pi * mm(trap_heater_od_mm)
    mcm_heater_length = mcm_heater_length_in * ureg['in']
    turns = (mcm_heater_length / heater_body_circumference).to('dimensionless'
        ).magnitude

    print(f'The heater cable will need to be wound {np.ceil(turns):.0f} times '
        'around the trap holder.'
    )

    # 1=all of the surface of the heater would be covered in adjacent turns.
    # >1 is not possible without overlap (which McMaster cautions against).
    winding_density = (heater_thickness * turns / mcm_heater_length).to(
        'dimensionless').magnitude

    # TODO is there any density short of something causing overlap that might
    # lead to failure? assuming not
    print('Winding density (close to 1 may mean overlap): '
        f'{winding_density:.3f}'
    )

    # TODO TODO need to check heater tube thickness likely enough (including at
    # maximum temperature!!) to withstand the desired pressure? or is even
    # 3mm Al in vast excess here?
    # The yield strength of most Al seems to drop to about 1/6th initial by
    # 350C, from "Overview of aluminum alloy mechanical properties during and
    # after fires", but "Barlow's formula" for expected bulge strength still
    # seems to be maybe ~1000psi using a ~0.5in OD and ~2.5mm wall.


# TODO TODO would it be sufficient to *just* raise the inlet gas up to the
# recommended desorption temperature, or would we also want to heat the
# trap to the same / an intermediate temperature? see patents?

# From all of the following product selections at fischersci.com (min - max):
# - Restek[tm] Instrument Grade Tubing
#   1.016 - 5.334 (ID in mm)
# - Restek[tm] Stainless-Steel HPLC Column Tubing
#   2.0828 - 4.5974
# - Restem[tm] Instrument-Grade Stainless Steel Tubing
#   0.254 - 5.334
# - Thermo Scientific[tm] 316 Stainless Steel Capillary Tubing for HPLC
#   in 5-Foot Coils
#   0.127 - 1.1684
# - Thermo Scientific[tm] 316 Stainless Steel Capillary Tubing for HPLC
#   0.127 - 0.254
# TODO make sure there are enough fittings available for all dimensions
# of the tubing i ultimately select
min_gas_heating_tube_id_mm = 0.127
max_gas_heating_tube_id_mm = 5.334

# TODO try nitrogen too? some reason it wouldn't work well for desorption?
# same kind of van deemter reasons as it's not used for GC carrier gas?
def min_gas_heater_power(gas='helium', max_vol_flow_l_per_m=2.0,
    pipe_diams_mm=None, pipe_lengths_mm=None, pipe_names=None):
    """
    max_vol_flow_l_per_m is the maximum volumetic flow desired at the
        output of the system (into ambient pressure).

    pipe_*_mm should be None or same-length iterables of those values,
        with values for upstream pipes listed first.
    """
    # TODO try diff values for this. need to be longer?
    gas_heating_tube_length_mm = 254

    if pipe_diams_mm is None:
        assert pipe_lengths_mm is None and pipe_names is None

        pipe_names = ['gas heating', 'sample trap', 'cooling']

        cooling_tube_length_mm = cooling_tube_length.to('mm').magnitude
        pipe_diams_mm = (
            gas_heating_tube_length_mm,
            trap_tube_length_mm,
            cooling_tube_length_mm
        )

    assert len(pipe_diams_mm) == len(pipe_lengths_mm)
    if pipe_names is not None:
        assert len(pipe_names) == len(pipe_diams_mm)

    import ipdb; ipdb.set_trace()
    
    # TODO start from exit pipe (assuming exiting to room pressure?
    # reasonable?) and calculate all pressure differentials up to mass flow
    # controller, to find pressure at exit of mass flow controller (or is there
    # a better way to approach this calculation?)

    # TODO maybe check reynolds in each pipe section (can they roughly
    # be considered independently?) to check laminar assumption OK?
    # (fail if not, unless supporting)

    # TODO incorporate assumption that desired (VOLUMETRIC) flow is
    # achieved at the exit of the system (mass flow everywhere in series
    # must be constant, but vol flow can vary as density changes with
    # pressure)

    #initial_pressure_psi = 

    import ipdb; ipdb.set_trace()


def main():
    #cooling()

    min_trap_heater_power()

    min_gas_heater_power()


if __name__ == '__main__':
    main()

