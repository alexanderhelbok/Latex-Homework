from astroquery.ned import Ned
from astropy.table import Table, vstack
import pandas as pd
import matplotlib.pyplot as plt

data = pd.DataFrame(columns=['Object Name', 'Velocity', 'Redshift', 'RA', 'DEC'])

# list of nebulas hubble used in his 1929 paper
objects = [221, 224, 598, 6822, 5457, 4736, 5194, 4449, 4214, 3031, 3627, 4826, 5236, 1068, 5055, 7331, 4258, 4151, 4382, 4472, 4486, 4649]

# Query NED for NGC 224 and NGC 225
for obj in objects:
    result_table = Ned.query_object('NGC ' + str(obj))
    # print(result_table['Object Name', 'Velocity', 'Redshift', 'RA', 'DEC'])
    astropy_df = result_table['Velocity', 'Redshift', 'RA', 'DEC'].to_pandas()
    data = pd.concat([data, astropy_df])
    # data = data.append(result_table['Velocity', 'Redshift', 'RA', 'DEC'], ignore_index=True)

data['Object Name'] = [f"NGC {Obj}" for Obj in objects]

print(data)
# extract velocity, redshift, and RA/DEC
# print(result_table)
# print(result_table['Velocity', 'Redshift', 'RA', 'DEC'])
# data = pd.DataFrame(result_table['Velocity', 'Redshift', 'RA', 'DEC'])


# plot velocity vs redshift
plt.scatter(data['Redshift'], data['Velocity'])
plt.ylabel('Velocity (km/s)')
plt.xlabel('Redshift')

plt.show()
