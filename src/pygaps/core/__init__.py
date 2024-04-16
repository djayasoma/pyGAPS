# pylint: disable=W0614,W0401,W0611,W0622
# flake8: noqa
import pygaps.parsing as pgp
from pygaps.units.converter_mode import c_isotherm_type


isotherm = pgp.isotherm_from_aif('./excess.aif')
new_isotherm = c_isotherm_type(
    iso = isotherm,
    total_pore_volume = 1.06,
    skeletal_density = 2.3,
    mode_from = 'excess',
    mode_to = 'total',
)
    
print(isotherm)
print(isotherm.data_raw)
print(new_isotherm.data_raw)
