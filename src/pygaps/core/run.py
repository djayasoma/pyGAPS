import pygaps.parsing as pgp
#import pygaps.parsing as pgg
import pygaps as pg


isotherm = pgp.isotherm_from_aif('./excess.aif')
#pgg.plot_iso(isotherm)

isotherm.material.skeletal_density = 2.0
isotherm.material.total_pore_volume = 1.3

new_isotherm = pg.PointIsotherm(
    pressure=isotherm.pressure(),
    loading=isotherm.loading(),
    temperature=isotherm.temperature,
    material=isotherm.material,
    adsorbate=isotherm.adsorbate,
    isotherm_type='excess',
    )
    
new_isotherm.c_isotherm_type(mode_to='total')
