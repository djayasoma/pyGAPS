"""Perform conversions between different variables used."""

from pygaps.units.converter_unit import _MASS_UNITS
from pygaps.units.converter_unit import _MOLAR_UNITS
from pygaps.units.converter_unit import _PRESSURE_UNITS
from pygaps.units.converter_unit import _TEMPERATURE_UNITS
from pygaps.units.converter_unit import _VOLUME_UNITS
from pygaps.units.converter_unit import _check_unit
from pygaps.units.converter_unit import c_unit
from pygaps.utilities.exceptions import ParameterError

import pygaps as pg
import pygaps.parsing as pgp
import pygaps.graphing as pgg
from pygaps.core.modelisotherm import ModelIsotherm
from pygaps.core.pointisotherm import PointIsotherm

from CoolProp.CoolProp import PropsSI

_PRESSURE_MODE = {
    "absolute": _PRESSURE_UNITS,
    "relative": None,
    "relative%": None,
}

_LOADING_MODE = {
    "mass": _MASS_UNITS,
    "volume_gas": _VOLUME_UNITS,
    "volume_liquid": _VOLUME_UNITS,
    "molar": _MOLAR_UNITS,
    "percent": None,
    "fraction": None,
}

_MATERIAL_MODE = {
    "mass": _MASS_UNITS,
    "volume": _VOLUME_UNITS,
    "molar": _MOLAR_UNITS,
}

_ISOTHERM_TYPE_MODE = {
    "total": None, 
    "excess": None,
    "net": None,
}


def _check_basis(basis, bases, btype):
    if not basis:
        raise ParameterError(
            "Specify the base to convert."
            f"Viable options are {list(bases.keys())}"
        )
    if basis not in bases:
        raise ParameterError(
            f"Basis selected for {btype} ({basis}) is not an option. "
            f"Viable options are {list(bases.keys())}"
        )

def c_pressure(
    value: float,
    mode_from: str,
    mode_to: str,
    unit_from: str,
    unit_to: str,
    adsorbate=None,
    temp: float = None,
):
    """
    Convert pressure units and modes.

    Adsorbate name and temperature have to be
    specified when converting between modes.

    Parameters
    ----------
    value : float
        The value to convert.
    mode_from : str
        Whether to convert from a mode.
    mode_to: str
        Whether to convert to a mode.
    unit_from : str
        Unit from which to convert.
    unit_to : str
        Unit to which to convert.
    adsorbate : Adsorbate, optional
        Adsorbate on which the pressure is to be
        converted. Required for mode change.
    temp : float, optional
        Temperature at which the pressure is measured, in K.
        Required for mode changes to relative pressure.

    Returns
    -------
    float
        Pressure converted as requested.

    Raises
    ------
    ``ParameterError``
        If the mode selected is not an option.
    """
    _check_basis(mode_from, _PRESSURE_MODE, 'pressure')
    _check_basis(mode_to, _PRESSURE_MODE, 'pressure')

    if mode_from != mode_to:

        unit = None
        sign = 1
        factor = 1

        # Now go through various global options
        if "absolute" in [mode_to, mode_from]:
            if mode_to == "absolute":
                _check_unit(unit_to, _PRESSURE_UNITS, 'pressure')
                unit = unit_to
                sign = 1

            if mode_from == "absolute":
                _check_unit(unit_from, _PRESSURE_UNITS, 'pressure')
                unit = unit_from
                sign = -1

            if not temp:
                raise ParameterError("A temperature is required for this conversion.")

            factor = adsorbate.saturation_pressure(temp, unit=unit)

            if "relative%" in [mode_to, mode_from]:
                factor = factor / 100

        elif mode_to in ["relative", "relative%"]:
            factor = 100
            if mode_to == "relative%":
                sign = 1
            elif mode_to == "relative":
                sign = -1

        return value * factor**sign

    # convert just units in absolute mode
    elif unit_to and mode_from == 'absolute':
        return c_unit(_PRESSURE_MODE[mode_from], value, unit_from, unit_to)

    # otherwise no change
    return value

def density(
    p: float,
    T: float,
    adsorbate: str,
):
    try:
        return PropsSI(
            'D',
            'T', T,
            'P|gas', p,
            adsorbate
        )
    except ValueError as e:
        print(
            f'Maybe one of your variables is the wrong type?\n'
            f'p:\t{type(p)}\n'
            f'T:\t{type(T)}\n'
            f'adsorbate:\t{type(adsorbate)}'
        )
        
def total_molar(
    density,
    excess_loading,
    total_pore_volume,
    molar_mass,
):
    excess_mass = excess_loading * molar_mass
    total_mass = excess_mass + (density * total_pore_volume)
    total_molar = total_mass / molar_mass
    return total_molar

def totexcess_molar(
    density,
    total_loading,
    total_pore_volume,
    molar_mass,
):
    total_mass = total_loading * molar_mass
    excess_mass = total_mass - (density * total_pore_volume)
    excess_molar = excess_mass / molar_mass
    return excess_molar

def net_molar(
    density,
    excess_loading,
    skeletal_density,
    molar_mass,
):
    excess_mass = excess_loading * molar_mass
    skeletal_volume = (1/skeletal_density)
    net_mass = excess_mass - (density * skeletal_volume)
    net_molar = net_mass / molar_mass
    return net_molar

def netexcess_molar(
    density,
    net_loading,
    skeletal_density,
    molar_mass,
):
    net_mass = net_loading * molar_mass
    skeletal_volume = (1/skeletal_density)
    excess_mass = net_mass + (density * skeletal_volume)
    excess_molar = excess_mass / molar_mass
    return excess_molar

def adsorption(
    isotherm: "ModelIsotherm|PointIsotherm",
    total_pore_volume: float,
    skeletal_density: float,
    mode_from: str,
    mode_to: str,
):

    isotherm.convert_loading(basis_to='molar', unit_to='mol')
    isotherm.convert_material(basis_to='mass', unit_to='kg')
    isotherm.convert_pressure(mode_to='absolute', unit_to='Pa')
    isotherm.convert_temperature(unit_to='K')

    adsorbate = isotherm.adsorbate
    molar_mass = isotherm.adsorbate.molar_mass()
    initial_loading = list(isotherm.loading())
    temperature = isotherm.temperature
    
    final_loading = []
    pressure = []
   
    for n in initial_loading:
        p = float(isotherm.pressure_at(n))
        if (n <= 0 or p<=0):
            continue

        d = density(
            p,
            temperature,
            str(adsorbate)
            )
        if mode_from == 'excess': #excess to total
            if mode_to == 'total':
                total = total_molar(
                    d, n,
                    total_pore_volume, molar_mass,)
                final_loading.append(total)
                print('total')
            elif mode_to == 'net':
                net = net_molar(
                    d, n,
                    skeletal_density, molar_mass,)
                final_loading.append(net)
                print('net')
        
        if mode_from == 'total': #excess to total
            if mode_to == 'excess':
                excess = totexcess_molar(
                    d, n,
                    total_pore_volume, molar_mass,)
                final_loading.append(excess)
                print('excess')
            elif mode_to == 'net':
                excess = totexcess_molar(
                    d, n,
                    total_pore_volume, molar_mass,)
                #return excess_molar
                excessnet = net_molar(
                    d, excess,
                    skeletal_density, molar_mass,)
                final_loading.append(excessnet)
                print('net')
        
        if mode_from == 'net': #excess to total
            if mode_to == 'total':
                excess = netexcess_molar(
                    d, n,
                    skeletal_density, molar_mass,)
                excesstot = total_molar(
                    d, excess,
                    total_pore_volume, molar_mass,)
                final_loading.append(excesstot)
                print('total')
            elif mode_to == 'excess':
                excess = netexcess_molar(
                    d, n,
                    skeletal_density, molar_mass,)
                final_loading.append(excess)
                print('excess')
        
        pressure.append(p)
            
    ads_isotherm = pg.PointIsotherm(
        pressure=pressure,
        loading=final_loading,
        
        material=isotherm.material,
        adsorbate=str(adsorbate),

        temperature=isotherm.temperature,
        temperature_unit='K',

        pressure_unit='Pa',
        pressure_mode='absolute',
        material_unit='kg',
        material_basis='mass',
        loading_unit='mol',
        loading_basis='molar',
    )
    return ads_isotherm


def c_isotherm_type(
    isotherm: "ModelIsotherm|PointIsotherm",
    total_pore_volume: float,
    skeletal_density: float,
    isotherm_type: str,
    mode_to: str,
):
    #import matplotlib.pyplot as plt
    _check_basis(isotherm_type, _ISOTHERM_TYPE_MODE, 'isotherm_type')
    _check_basis(mode_to, _ISOTHERM_TYPE_MODE, 'isotherm_type')

    if isotherm_type != mode_to:
        #isotherm = pgp.isotherm_from_aif('./excess.aif')
        ads_isotherm = adsorption(isotherm, total_pore_volume, skeletal_density, isotherm_type, mode_to)
        for iso in [isotherm, ads_isotherm]:
            iso.convert(
                pressure_unit='bar',
                loading_unit='mmol',
                material_unit='g',
            )
        #pgp.isotherm_to_aif(ads_isotherm,'result2.aif')
                  
        #pgg.plot_iso(
        #    [isotherm, ads_isotherm]
        #    )
        #plt.legend(['1st', '2nd'])
        #plt.show()
        #print(iso)
        #print(isotherm.data_raw)

def c_loading(
    value: float,
    basis_from: str,
    basis_to: str,
    unit_from: str,
    unit_to: str,
    adsorbate=None,
    temp: float = None,
    basis_material: str = None,
    unit_material: str = None,
):
    """
    Convert loading units and basis.

    Adsorbate name and temperature have to be
    specified when converting between basis.

    Parameters
    ----------
    value : float
        The value to convert.
    basis_from : str
        Whether to convert from a basis.
    basis_to: str
        Whether to convert to a basis.
    unit_from : str
        Unit from which to convert.
    unit_to : str
        Unit to which to convert.
    adsorbate : str, optional
        Adsorbate for which the pressure is to be converted.
        Only required for some conversions.
    temp : float, optional
        Temperature at which the loading is measured, in K.
        Only required for some conversions.
    basis_material : str, optional
        The basis of the material.
        Only required for conversions involving percentage/fraction.
    unit_material : str, optional
        The unit of the material.
        Only required for conversions involving percentage/fraction.

    Returns
    -------
    float
        Loading converted as requested.

    Raises
    ------
    ``ParameterError``
        If the mode selected is not an option.
    """
    _check_basis(basis_from, _LOADING_MODE, 'loading')
    _check_basis(basis_to, _LOADING_MODE, 'loading')

    if basis_from != basis_to:

        if _LOADING_MODE[basis_to]:
            _check_unit(unit_to, _LOADING_MODE[basis_to], 'loading')

        if _LOADING_MODE[basis_from]:
            _check_unit(unit_from, _LOADING_MODE[basis_from], 'loading')

        _basis_from = basis_from
        _unit_from = unit_from
        _basis_to = basis_to
        _unit_to = unit_to

        constant = 1
        sign = 1
        factor = 1

        bf_b = basis_from in ['percent', 'fraction']
        bt_b = basis_to in ['percent', 'fraction']
        if bf_b or bt_b:
            # if both are percent/fraction we do not need conversion
            if bf_b and bt_b:
                # we know they are different, so one must be percent and one fraction
                if basis_from == 'percent':
                    return value / 100
                else:
                    return value * 100
            else:
                # convert from physical -> percent/fraction
                if basis_material == 'volume':
                    basis_material = 'volume_liquid'
                if bf_b:
                    _basis_from = basis_material
                    _unit_from = unit_material
                    if basis_from == 'percent':
                        factor = 0.01
                elif bt_b:
                    _basis_to = basis_material
                    _unit_to = unit_material
                    if basis_to == 'percent':
                        factor = 100

        # TODO: this list grew, should find a better way
        if _basis_from == 'mass':
            if _basis_to == 'volume_gas':
                constant = adsorbate.gas_density(temp=temp)
                sign = -1
            elif _basis_to == 'volume_liquid':
                constant = adsorbate.liquid_density(temp=temp)
                sign = -1
            elif _basis_to == 'molar':
                constant = adsorbate.molar_mass()
                sign = -1
        elif _basis_from == 'volume_gas':
            if _basis_to == 'mass':
                constant = adsorbate.gas_density(temp=temp)
                sign = 1
            elif _basis_to == 'molar':
                constant = adsorbate.gas_molar_density(temp=temp)
                sign = 1
            elif _basis_to == 'volume_liquid':
                constant = adsorbate.gas_molar_density(temp=temp) /\
                    adsorbate.liquid_molar_density(temp=temp)
                sign = 1
        elif _basis_from == 'volume_liquid':
            if _basis_to == 'mass':
                constant = adsorbate.liquid_density(temp=temp)
                sign = 1
            elif _basis_to == 'molar':
                constant = adsorbate.liquid_molar_density(temp=temp)
                sign = 1
            elif _basis_to == 'volume_gas':
                constant = adsorbate.gas_molar_density(temp=temp) /\
                    adsorbate.liquid_molar_density(temp=temp)
                sign = -1
        elif _basis_from == 'molar':
            if _basis_to == 'mass':
                constant = adsorbate.molar_mass()
                sign = 1
            elif _basis_to == 'volume_gas':
                constant = adsorbate.gas_molar_density(temp=temp)
                sign = -1
            elif _basis_to == 'volume_liquid':
                constant = adsorbate.liquid_molar_density(temp=temp)
                sign = -1

        return (
            value * _LOADING_MODE[_basis_from][_unit_from] * factor * constant**sign /
            _LOADING_MODE[_basis_to][_unit_to]
        )

    if unit_to and unit_from != unit_to:
        return c_unit(_LOADING_MODE[basis_from], value, unit_from, unit_to)

    return value


def c_material(
    value: float,
    basis_from: str,
    basis_to: str,
    unit_from: str,
    unit_to: str,
    material=None,
):
    """
    Convert material units and basis.

    The name of the material has to be
    specified when converting between basis.

    Parameters
    ----------
    value : float
        The value to convert.
    basis_from : str
        Whether to convert from a basis.
    basis_to: str
        Whether to convert to a basis.
    unit_from : str
        Unit from which to convert.
    unit_to : str
        Unit to which to convert.
    material : str
        Name of the material on which the value is based.

    Returns
    -------
    float
        Loading converted as requested.

    Raises
    ------
    ``ParameterError``
        If the mode selected is not an option.

    """
    _check_basis(basis_from, _MATERIAL_MODE, 'material')
    _check_basis(basis_to, _MATERIAL_MODE, 'material')

    if basis_from != basis_to:

        if (basis_from in ['percent', 'fraction'] or basis_to in ['percent', 'fraction']):
            raise ParameterError(
                "If you want to convert to/from fraction/percent,"
                " convert using loading, not adsorbate."
            )

        _check_unit(unit_to, _MATERIAL_MODE[basis_to], 'material')
        _check_unit(unit_from, _MATERIAL_MODE[basis_from], 'material')

        constant = 1
        sign = 1

        if basis_from == 'mass':
            if basis_to == 'volume':
                constant = material.density
                sign = -1
            elif basis_to == 'molar':
                constant = material.molar_mass
                sign = -1
        elif basis_from == 'volume':
            if basis_to == 'mass':
                constant = material.density
                sign = 1
            elif basis_to == 'molar':
                constant = material.density / material.molar_mass
                sign = 1
        elif basis_from == 'molar':
            if basis_to == 'mass':
                constant = material.molar_mass
                sign = 1
            elif basis_to == 'volume':
                constant = material.density / material.molar_mass
                sign = -1

        return (
            value / _MATERIAL_MODE[basis_from][unit_from] / constant**sign *
            _MATERIAL_MODE[basis_to][unit_to]
        )

    if unit_to and unit_from != unit_to:
        return c_unit(_MATERIAL_MODE[basis_from], value, unit_from, unit_to, sign=-1)

    return value


def c_temperature(
    value: float,
    unit_from: str,
    unit_to: str,
):
    """
    Convert temperatures.

    Parameters
    ----------
    value : float
        The value to convert.
    unit_from : str
        Unit from which to convert.
    unit_to : str
        Unit to which to convert.

    Returns
    -------
    float
        Temperature converted as requested.

    Raises
    ------
    ``ParameterError``
        If the unit selected is not an option.

    """
    if unit_to and "c" in unit_to.lower():
        unit_to = "°C"
    if unit_from and "c" in unit_from.lower():
        unit_from = "°C"
    _check_unit(unit_to, _TEMPERATURE_UNITS, 'temperature')
    _check_unit(unit_from, _TEMPERATURE_UNITS, 'temperature')

    if unit_from == unit_to:
        return value

    return value - _TEMPERATURE_UNITS[unit_to]
