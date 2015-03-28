
from gammapy.datasets import fetch_fermi_catalog, fetch_fermi_extended_sources

#get point source catalog
#fetch_fermi_catalog('2FGL')
fetch_fermi_catalog('3FGL')

#get extended source catalog
#fetch_fermi_extended_sources('2FGL')
fetch_fermi_extended_sources('3FGL')

#check 3FGL properties
n_hdu = len(fetch_fermi_catalog('3FGL'))
n_sources_point = len(fetch_fermi_catalog('3FGL', 'LAT_Point_Source_Catalog'))
n_sources_ext = len(fetch_fermi_extended_sources('3FGL'))

print "3FGL catalog: %i HDUs, %i point sources, %i extended sources" %(n_hdu, n_sources_point, n_sources_ext)
