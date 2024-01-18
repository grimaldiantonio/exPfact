"""
Takes as input a cluster file from DnX and generates an assignment file (coverage map)
and a .dexp file.

*** Back-exchange correction is with respect to columns MaxUptake
*** Different charge states are averaged together
"""

import pandas as pd
import numpy as np

infile = "CTB_cluster.csv"
state = "CTB+GM1os"

data = pd.read_csv(infile)

data['Dt'] = data['Center'] * data['z'] - data['z']
data['CorrDt'] = 0

# Find coverage map - peptides and sequences
peptides = []
sequences = []
for i in range(len(data)):
    if (data['Start'][i], data['End'][i]) not in peptides:
        peptides.append((data['Start'][i], data['End'][i]))
        if data['Sequence'][i] not in sequences:
            sequences.append(data['Sequence'][i])

# Write assignment file
with open(infile.replace(".csv","_%s.ass" % state), 'w') as f:
    for i in range(len(peptides)):
        f.write("%d %d %d %s\n" % (i+1, peptides[i][0], peptides[i][1], sequences[i]))

times = np.sort(list(set(data['Exposure'])))

dexp = np.zeros((len(peptides)+1, len(times)))
dexp[0] = times
for pep in range(1, len(peptides)+1):
    subdata = data[(data['State'] == state) &\
                   (data['Start'] == peptides[pep-1][0]) &\
                   (data['End'] == peptides[pep-1][1])]        
    control = np.average(subdata[subdata['Exposure']==0]['Dt'])
    maxUptake = list(set(subdata['MaxUptake']))[0]
    for j in range(len(times)):
        data_t = subdata[subdata['Exposure'] == times[j]]
        Dt = (data_t['Dt'] - control) / maxUptake
        uptake = np.average(Dt)
        dexp[pep][j] = uptake
np.savetxt(infile.replace(".csv", "_%s.dexp" % state), np.transpose(dexp), fmt = "%5.5f")
