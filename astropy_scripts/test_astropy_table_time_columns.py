from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3
import numpy as np
from astropy.table import Table, Column
from astropy.time import Time, TimeDelta
from astropy.units import Quantity

# create an astrotable with a column of integers
n_obs = 10 # number of rows
n_obs_start = 1 # start counting at 1
astro_table = Table()
col_obs_id = Column(name='obs_id', data=np.arange(n_obs_start, n_obs_start + n_obs))
astro_table.add_column(col_obs_id)

# try Time and TimeDelta columns

# start time
# random points between the start of 2010 and the end of 2014
datestart = Time('2010-01-01 00:00:00', format='iso', scale='utc')
dateend = Time('2015-01-01 00:00:00', format='iso', scale='utc')
time_start = Time((dateend.mjd - datestart.mjd)*np.random.random(len(col_obs_id)) + datestart.mjd, format='mjd', scale='utc').iso
col_time_start = Column(name='TSTART', data=time_start)
astro_table.add_column(col_time_start)

# on time: 30 min
ontime = TimeDelta(30.*60.*np.ones_like(col_obs_id.data), format='sec')
#col_ontime = Column(name='ONTIME', data=ontime) #faulty!
#astro_table.add_column(col_ontime)
astro_table['ONTIME'] = ontime #works fine, but still no unit showing in the column
# ref: https://github.com/astropy/astropy/issues/3832

# try Quantity instead of TimeDelta quantity

# livetime: 25 min
livetime = Quantity(25.*60.*np.ones_like(col_obs_id.data), 'second')
col_livetime = Column(name='LIVETIME', data=livetime)
astro_table.add_column(col_livetime)

print("astro_table")
print(astro_table)

print()

print("type time col", type(astro_table['TSTART']))
print("type timedelta col", type(astro_table['ONTIME']))
print("type quantity col", type(astro_table['LIVETIME']))

# read columns from the table
#import IPython; IPython.embed() # before error line
time_col = Time(astro_table['TSTART']) # works fine
timedelta_col = TimeDelta(astro_table['ONTIME']) # gives an error!
quantity_col = TimeDelta(astro_table['LIVETIME']) # works fine

print()

# better output if running ipython!
print(time_col.__dict__)
print(timedelta_col.__dict__) # faulty!
print(quantity_col.__dict__)
