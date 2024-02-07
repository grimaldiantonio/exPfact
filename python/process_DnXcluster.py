"""
Takes as input a cluster file from DnX and generates an assignment file (coverage map)
and a .dexp file.

*** Back-exchange correction is with respect to columns MaxUptake
*** Different charge states are averaged together
"""

import pandas as pd
import numpy as np
import sys

if __name__ == '__main__':
    
    if len(sys.argv) > 1:
        infile = sys.argv[1]
    else:
        print("Error\nUsage: process_DnXcluster.py <filename.csv>")
        exit()

    alldata = pd.read_csv(infile)
    print("Converting DynamX cluster %s file into exPfact input files\n" % infile)

    alldata['Dt'] = alldata['Center'] * alldata['z'] - alldata['z']
    alldata["Time"] = round(alldata["Exposure"], 2)
    alldata['CorrDt'] = 0

    states = list(set(alldata["State"]))
    print("\nThe following experimental conditions have been identified:")
    for i in range(len(states)):
        print("  %d) %s" % (i+1, states[i]))
    
    times = np.sort(list(set(alldata['Exposure'])))
    print("\nThe following labelling times have been identified:")
    for i in range(len(times)):
        print("  %d) %5.1f min" % (i+1, times[i]))
        
    print("\nExPfact needs deuterium uptake values to be normalized between 0 and 1.")
    print("Two normalization methods are available:")
    print("  1) Using the theoretical maximum uptake (column MaxUptake in the DynamX cluster file, multiplied by D2O percentage)")
    print("  2) Using the fully deuterated sample (the longer exposure time is used as fully deuterated sample)")
    print("Please select the method that you prefer [1/2]:")
    norm_method = input()

    if norm_method == "1":
        print("\nPlease state the D2O percentage [0 - 100]:")
        d2o_perc = int(input())
    elif norm_method == "2":
        back_exs = []
        allpeptides = list(set(alldata["Sequence"]))
        for peptide in allpeptides:
            subdata = alldata[(alldata["Sequence"] == peptide) &\
                              (alldata["Exposure"] == times[-1])].reset_index()
            fddata = list(subdata["Dt"])
            if len(fddata) > 1:
                mhp = np.average(list(subdata["MHP"]))
                maxu = np.average(list(subdata["MaxUptake"]))
                fd = np.average(fddata)
                back_ex = (fd - mhp) / maxu 
                back_exs.append(back_ex)
        back_ex_perc = np.average(back_exs) * 100
        print("\nAverage back-exchange plateau: %5.2f" % back_ex_perc)
        print("calculated as ratio between maximally labelled sample and theoretical maximal uptake")
        print("when not available, fully deuterated sample estimated using the average back exchange percentage")
    else:
        print("\nNormalization method not valid")
        exit()
    print("\n")

    for state in states:     
        data = alldata[alldata["State"] == state].reset_index()
        
        # Find coverage map - peptides and sequences
        peptides = []
        sequences = []
        for i in range(len(data)):
            if (data['Start'][i], data['End'][i]) not in peptides:
                peptides.append((data['Start'][i], data['End'][i]))
                if data['Sequence'][i] not in sequences:
                    sequences.append(data['Sequence'][i])
        
        print("State: %s" % state)
        print("  Number of peptides: %d" % len(peptides))
        
        # Write assignment file
        with open(infile.replace(".csv","_%s.list" % state), 'w') as f:
            for i in range(len(peptides)):
                f.write("%d %d %d %s\n" % (i+1, peptides[i][0], peptides[i][1], sequences[i]))
        
        times = np.sort(list(set(data['Time'])))
        
        dexp = np.zeros((len(peptides)+1, len(times)))
        dexp[0] = times / 60
        for pep in range(1, len(peptides)+1):
            subdata = data[(data['Start'] == peptides[pep-1][0]) &\
                           (data['End'] == peptides[pep-1][1])]        
            control = np.average(subdata[subdata['Time'] == times[0]]['Dt'])
            
            if norm_method == "1":
                maxUptake = d2o_perc/100 * list(set(subdata['MaxUptake']))[0]
            elif norm_method == "2":
                fd_control = subdata[subdata['Time'] == times[-1]]['Dt']
                if len(fd_control) > 0:
                    maxUptake = np.average(fd_control)
                else:
                    maxUptake = back_ex_perc/100 * list(set(subdata['MaxUptake']))[0]
            else:
                print("Normalization method not valid!")
                exit()
    
            for j in range(len(times)):
                data_t = subdata[subdata['Time'] == times[j]]
                if norm_method == "1":
                    Dt = (np.average(data_t['Dt']) - control) / maxUptake
                elif norm_method == "2":
                    Dt = (data_t['Dt'] - control) / (maxUptake - control)
                if len(Dt) > 0:
                    uptake = np.average(Dt)
                else:
                    print("  Note: missing time point %5.1f for peptide %s in state %s" % (times[j], sequences[pep-1], state))
                    uptake = "NaN"
                
                if uptake < 0:
                    dexp[pep][j] = 0
                elif uptake > 1:
                    dexp[pep][j] = 1
                else:
                    dexp[pep][j] = uptake
                    
        if norm_method == "2":
            np.savetxt(infile.replace(".csv", "_%s.dexp" % state), np.transpose(dexp)[:-1], fmt = "%5.5f")
        else:
            np.savetxt(infile.replace(".csv", "_%s.dexp" % state), np.transpose(dexp), fmt = "%5.5f")
