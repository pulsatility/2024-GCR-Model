#%% Loading libraries
import os
import math
import random
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from collections import Counter
import seaborn as sns
from statistics import mean, median, stdev
from scipy import stats


#%% Define simulation dataset directory
"""
Define simulation dataset directory
"""

# Main result: Dataset_V50_33
# Fig. S1A-S1B: Dataset_V50_31
# Fig. S1C-S1D: Dataset_V50_32
# Fig. S4: Dataset_V51_01
# Fig. S5: Dataset_V50_02
 
# Manually set the current directory to where the saved .txt files are stored, then obtain the directory string. For instance: C:\Users\qzhan31\GCR_saved_files\Dataset_V50_33
cellData_dir = os.getcwd()


#%% Define variables and functions
"""
Define variables and functions
"""

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
font = fm.FontProperties(family='Arial')
save_interval = 15
xOffset = 25

# DZ, midline, and LZ boundaries
DZxCOM = 100 + xOffset
midxCOM = 100 + xOffset
LZxCOM = 100 + xOffset

variables_not_to_display = ['affinityThreshold', 'cellID', 'motherID', 'seederID', 'numMutation', \
                          'deathMutation', 'Alive', 'selectionDecision', 'Cytokinesis', 'plasmaCell', \
                              'thresholdVolumeDivision', 'targetVolume']
  


# Function that evaluates whether the variable is a float
def is_float(element):
    try:
        float(element)
        return True
    except ValueError:
        return False

# Function that reads in all lines in a .txt file
def readFile(id):
    fileName = id + ".txt"
    with open(fileName) as f:
        lines = f.readlines()
    return lines


#%% Simulation data pre-processing
"""
Simulation data pre-processing
"""

# Obtain names of all .txt files in the specified folder and store them as string list and string tuple
fileList = []
for fileName in os.listdir(cellData_dir):
    if fileName[-4:] == '.txt':
        fileList.append(fileName[:-4])
fileTuple = tuple(fileList) #Convert to tuple as tuple is read-only.

# Find total number of files (B cells)
numFiles = len(fileTuple)
print("Total number of files: " + str(numFiles))

# Extract header using the first line in the first file
header = readFile(fileTuple[0])[0].split()

# Obtain the index of each variable name in header
for variable_name in header:
    globals()['idx_'+variable_name] = header.index(variable_name)
    

# Merge time-course data of all cells vertically into a list and a tuple
data_all_list = []
for file in fileList:
    for line in readFile(file):  # Process one line at a time in current file
        lineList = line.split()  # convert a line to a list of strings
        if is_float(lineList[idx_mcs]):  # check that the line is not a header (which is supposed to exist only in seeder cell files) by using mcs
            #data_all.append([lineList[0]] + lineList[idx_initialAffinity:-1]) #Use this if need to skip the AntibodySequence column in data. The double ]] is needed to convert lineList[0] back to a list.
            data_all_list.append(lineList)
data_all_tuple = tuple(data_all_list)

# Obtain all mcs points (including those when a cell dies, divides, differentiates) and unique, sorted mcs points 
mcs_all = [int(a[idx_mcs]) for a in data_all_tuple] #Extract all mcs points
mcs_set = list(set(mcs_all))
mcs_set.sort()
mcs_set_sorted = mcs_set # Unique, sorted list of mcs points
mcs_max = max(mcs_set_sorted)
mcs_regular_saves = list(range(0, mcs_max+1, save_interval))  # range() is exclusive
mcs_regular_save_last = mcs_regular_saves[-1]


# Find all seeder B cells
cellIDs_seeder = [i for i in fileTuple if i.startswith('0_')]


#%% Fig 2E, 9, S3

# Function to retrieve all lineage branches
def retrieveLineageBranches(currID, data, txt):
    for line in readFile(currID):  # read in current file
        data += [line.split()]
    txt += [currID]
    children = []
    exit = True
    for file in fileList_copy:  # find children
        if file.split("_")[1] == currID.split("_")[2]:
            children.append(file)
            exit = False
    for child in children:  # recurse on children
        if child in fileList_copy:
            fileList_copy.remove(child)
        retrieveLineageBranches(child, list(data), list(txt))
    if exit == True:  # at root
        lineages.append(txt)


# Retrieve all lineage branches for all seeder cells
fileList_copy = fileList[:]
if input("Need lineages? (Y/N): ") == 'Y':
    lineages = []
    for seeder in cellIDs_seeder:
        if seeder in fileList_copy:
            fileList_copy.remove(seeder) # Prevents roots from calling themselves as daughter cells because their parent id = their cell id
        retrieveLineageBranches(seeder, [], [])
    for i, lineage in enumerate(lineages):
        print(i)
        #print(lineage)


# Find lineage branches that contain a particular cellID
cellID = "36_14597_15414" #“25_6946_7422”， #"20_4023_4311"， #"36_14597_15414"， #"27_5253_5614"， #"21_3961_4127"
indice_chosenBranches = [lineages.index(branch) for branch in lineages if cellID in branch]
print(indice_chosenBranches)


# Choose one lineage branch and vertically combine data of all cells in the chosenBranch without the header
chosenBranch = lineages[14073] #14073 #25311 #38625 for Fig. 2E. #35030 for Fig. 9, #25621 #51135 for Fig. S3.  #5580 #9173 #15700
#chosenBranch = lineages[indice_chosenBranches[0]]
data_chosenBranch = []
for file in chosenBranch:
    for line in readFile(file):  # read in current file
        lineList = line.split()
        if is_float(lineList[idx_mcs]):
                data_chosenBranch.append(lineList)


"""
Generating Fig. 9, S3 - Single lineage branch trajectory
"""
# Automatically extract each variable's values into a list and plot time-course (Fig. 9 and S3)
for variable_name in header:
    if variable_name not in variables_not_to_display:
        idx = header.index(variable_name)
        globals()[variable_name] = [float(row[idx]) for row in data_chosenBranch]
        fig, ax = plt.subplots()
        plt.plot(np.divide(mcs[0:-1:2], 100), globals()[variable_name][0:-1:2], label=variable_name)
        plt.ylim(top=max(globals()[variable_name][0:-1:2])*1.15)
        plt.xlabel('Time (x100 mcs)',fontsize=16)
        #plt.ylabel(variable_name,fontsize=16)
        plt.title(variable_name,fontsize=18)
        plt.xticks(range(0, round(max(np.divide(mcs[0:-1:2], 100)))+10, 50))
        #plt.legend(fontsize=24, frameon=False, handlelength=0, loc='upper center', bbox_to_anchor=(0.47, 1.05))
        ax.tick_params(axis='both', direction='in', labelsize=20)
        plt.savefig('Fig. 10 - ' + variable_name + '.png', dpi=600)
        plt.show()
        

# To customize a particular plot
fig, ax = plt.subplots()
plt.plot(np.divide(mcs[0:-1:2], 100), xCOM[0:-1:2], label='xCOM')
plt.plot(np.divide(mcs[0:-1:2], 100), actualVolume[0:-1:2], label='Cell volume')
plt.xlabel('Time (x100 mcs)')
plt.legend()
plt.xticks(fontsize=24)
plt.xticks(range(0, round(max(np.divide(mcs[0:-1:2], 100)))+10, 50))
plt.yticks(fontsize=24)
plt.legend(fontsize=24, frameon=False, handlelength=0, loc='upper center', bbox_to_anchor=(0.47, 1.05))
plt.show()


"""
Generating Fig 2E - xCOM-yCOM scatter with mcs color
"""
hsv_modified = cm.get_cmap('hsv', 256) # modified hsv in 256 color class
newcmp = ListedColormap(hsv_modified(np.linspace(0.7, 0, 256))) #create new hsv colormaps in range of 0 (red) to 0.7 (blue)

fig, ax = plt.subplots()
plt.scatter(xCOM[0:-1:1], yCOM[0:-1:1], c=np.divide(mcs[0:-1:1], 100), s= 6,cmap=newcmp) #cmap='rainbow_r'
cbar = plt.colorbar(ticks=[0, 100, 200, 300], location='top', shrink=0.65)
for t in cbar.ax.get_xticklabels():
     t.set_fontsize(14)
ax.set_aspect(1.35)
plt.xlabel('xCOM',fontsize=16)
plt.ylabel('yCOM',fontsize=16)
plt.xticks(fontsize=16)
plt.xticks(range(0, 251, 50))
plt.yticks(fontsize=16)
plt.yticks(range(25, 200, 50))
plt.savefig('Fig. 2E.png', dpi=600)
plt.show()



#%% Population dynamics
"""
Obtain cell IDs of different fates
"""

# Obtain lists of cells that have died ("Alive" in last line is 0 in .txt file) or never died ("Alive" is always 1 in .txt file)
cellIDs_alive = []
cellIDs_alive_cytokinesis = []
cellIDs_alive_plasmaCell = []
cellIDs_alive_simulationEnd = []
cellIDs_alive_volumeVanish = []

cellIDs_dead = []
cellIDs_dead_deathMutation = []
cellIDs_dead_deathTimer = []
cellIDs_dead_deathTimer_neg_selected = []
cellIDs_dead_deathTimer_noTtouch = []
cellIDs_dead_deathTimer_noTtouch_DZ = []
cellIDs_dead_deathTimer_noTtouch_LZ = []

cellIDs_ready_to_move_to_LZ = []

for file in fileTuple:
    
    lines = readFile(file)  # read in all lines in current file
    firstlineList = lines[0].split()
    lastlineList = lines[-1].split()
    
    if is_float(lastlineList[idx_Alive]):
        if int(lastlineList[idx_Alive]) == 0:
            cellIDs_dead.append(file)
            if int(lastlineList[idx_deathMutation]) == 1:
                cellIDs_dead_deathMutation.append(file)
            else:
                cellIDs_dead_deathTimer.append(file)
                if int(lastlineList[idx_selectionDecision]) == 1:
                    cellIDs_dead_deathTimer_neg_selected.append(file)
                else:
                    cellIDs_dead_deathTimer_noTtouch.append(file)
                    if float(lastlineList[idx_xCOM]) <= midxCOM:
                        cellIDs_dead_deathTimer_noTtouch_DZ.append(file)
                    else:
                        cellIDs_dead_deathTimer_noTtouch_LZ.append(file)        
                
        elif int(lastlineList[idx_Alive]) == 1:
            cellIDs_alive.append(file)
            if int(lastlineList[idx_Cytokinesis]):
                cellIDs_alive_cytokinesis.append(file)
            elif int(lastlineList[idx_plasmaCell]) == 1:
                cellIDs_alive_plasmaCell.append(file)
            elif int(lastlineList[idx_mcs]) == mcs_max or int(lastlineList[idx_mcs]) == mcs_regular_save_last: #Some volumeVanishing cells may still terminate here
                cellIDs_alive_simulationEnd.append(file)
            else:
                cellIDs_alive_volumeVanish.append(file)
                
    if is_float(firstlineList[idx_cellCycleCommitted]): #This condition will exclude seeder cells which started in the midZone
        if float(firstlineList[idx_cellCycleCommitted]) == 0: #If a cell start with cellCycleCommitted = 0, it must be a seeder cell or cell that has existed the cell cylce ready to move to LZ if in DZ. These will include some that has exited the cell cycle but will die due to lethal mutation.
            cellIDs_ready_to_move_to_LZ.append(file)

# Check if all cell numbers add up
if len(cellIDs_alive)+len(cellIDs_dead) == numFiles and len(cellIDs_dead_deathMutation)+len(cellIDs_dead_deathTimer) == len(cellIDs_dead):
    print('All numbers add up')
else:
    print("Some numbers don't add up!!!")



#%% Fig. 3A, S4A, S5A - cell counts
"""
Fig. 3A, S4A, S5A - cell counts
"""
# Calculate number of existing cells in whole GC at regular save interval
mcs_all_regular = [i for i in mcs_all if i in mcs_regular_saves]
counts_all_at_mcs_regular = Counter(mcs_all_regular) #Obtain unique mcs and corresponding count
counts_all_at_mcs_regular = sorted(counts_all_at_mcs_regular.items()) #Sort by the key (mcs) and convert to a list of tuple
x_all = [tple[0] for tple in counts_all_at_mcs_regular]
y_all = [tple[1] for tple in counts_all_at_mcs_regular]

# Calculate number of existing cells in DZ at regular save interval
mcs_all_DZ = [int(a[idx_mcs]) for a in data_all_tuple if float(a[idx_xCOM])<DZxCOM ] #Extract all mcs points for DZ cells
mcs_all_DZ_regular = [i for i in mcs_all_DZ if i in mcs_regular_saves]
counts_all_DZ_at_mcs_regular = Counter(mcs_all_DZ_regular) #Obtain unique mcs and corresponding count
counts_all_DZ_at_mcs_regular = sorted(counts_all_DZ_at_mcs_regular.items()) #Sort by the key (mcs) and convert to a list of tuple
x_DZ = [tple[0] for tple in counts_all_DZ_at_mcs_regular]
y_DZ = [tple[1] for tple in counts_all_DZ_at_mcs_regular]

# Calculate number of existing cells in midZ at regular save interval
mcs_all_midZ = [int(a[idx_mcs]) for a in data_all_tuple if float(a[idx_xCOM])>=DZxCOM and float(a[idx_xCOM])<=LZxCOM ] #Extract all mcs points for DZ cells
mcs_all_midZ_regular = [i for i in mcs_all_midZ if i in mcs_regular_saves]
counts_all_midZ_at_mcs_regular = Counter(mcs_all_midZ_regular) #Obtain unique mcs and corresponding count
counts_all_midZ_at_mcs_regular = sorted(counts_all_midZ_at_mcs_regular.items()) #Sort by the key (mcs) and convert to a list of tuple
x_midZ = [tple[0] for tple in counts_all_midZ_at_mcs_regular]
y_midZ = [tple[1] for tple in counts_all_midZ_at_mcs_regular]

# Calculate number of existing cells in LZ at regular save interval
mcs_all_LZ = [int(a[idx_mcs]) for a in data_all_tuple if float(a[idx_xCOM])>LZxCOM ] #Extract all mcs points for DZ cells
mcs_all_LZ_regular = [i for i in mcs_all_LZ if i in mcs_regular_saves]
counts_all_LZ_at_mcs_regular = Counter(mcs_all_LZ_regular) #Obtain unique mcs and corresponding count
counts_all_LZ_at_mcs_regular = sorted(counts_all_LZ_at_mcs_regular.items()) #Sort by the key (mcs) and convert to a list of tuple
x_LZ = [tple[0] for tple in counts_all_LZ_at_mcs_regular]
y_LZ = [tple[1] for tple in counts_all_LZ_at_mcs_regular]

# Calculate DZ:LZ ratio. This method skips all mcs moments when either DZ or LZ has no B cells
DZ_LZ_ratio = []
mcs_for_ratio = []
for xDZ in x_DZ:
    if xDZ in x_LZ:
        DZ_LZ_ratio.append(y_DZ[x_DZ.index(xDZ)] / y_LZ[x_LZ.index(xDZ)])
        mcs_for_ratio.append(xDZ)
# DZ_LZ_ratio = [DZ/LZ for DZ, LZ in zip(y_DZ, y_LZ)] # Old code, only works if both have the same length

# Plot cell numbers vs mcs at regular saving intervals
fig, axL = plt.subplots()
axR = axL.twinx()  

axR.plot(np.divide(mcs_for_ratio,100), DZ_LZ_ratio, c="black", ls=":", label="DZ:LZ ratio")
axR.set_yscale("log")
axR.set_ylabel("DZ:LZ ratio", fontsize=16)

axL.plot(np.divide(x_all,100), np.divide(y_all,1000), c="blue", label="Total")
axL.plot(np.divide(x_DZ,100), np.divide(y_DZ,1000), c="red", label="DZ")
axL.plot(np.divide(x_LZ,100), np.divide(y_LZ,1000), c="green", label="LZ")
#axL.plot(x_midZ, y_midZ, c="gray", label="midZ")
axL.set_ylabel("Number of B cells", fontsize=16)
axL.tick_params(direction='in', labelsize=20)

hR, lR = axR.get_legend_handles_labels()
hL, lL = axL.get_legend_handles_labels()
axL.legend(hR+hL, lR+lL, loc='lower left', bbox_to_anchor=(0.052, 0.16))

axR.set_xlabel("Time (x100 mcs)", fontsize=16)
axR.set_title("Number of B cells",fontsize=16)
axR.tick_params(direction='in', labelsize=20)
plt.xticks(range(0,721,240))
plt.savefig('Fig. 3A.png', dpi=600, bbox_inches = "tight") # Change file name to Fig. S4A or Fig. S5A as needed
plt.show()



#%% Fig. 3B, 3C, 3E, 3F - Cell cycle committed and cell birth
"""
Fig. 3B, 3C - Cell cycle committed
"""

# Calculate and plot number of cells committed to cell cycle at regular save interval
#All
mcs_all_cellCycleCommitted = [int(a[idx_mcs]) for a in data_all_tuple if float(a[idx_cellCycleCommitted])==1 ] #Extract all mcs points when cellCycleCommitted=1
mcs_all_cellCycleCommitted_regular = [i for i in mcs_all_cellCycleCommitted if i in mcs_regular_saves]
counts_all_cellCycleCommitted_at_mcs_regular = Counter(mcs_all_cellCycleCommitted_regular)
counts_all_cellCycleCommitted_at_mcs_regular = sorted(counts_all_cellCycleCommitted_at_mcs_regular.items()) #Sort by the key (mcs) and convert to a list of tuple
x_all_cellCycleCommitted = [tple[0] for tple in counts_all_cellCycleCommitted_at_mcs_regular]
y_all_cellCycleCommitted = [tple[1] for tple in counts_all_cellCycleCommitted_at_mcs_regular]
y_all_allcells = [y_all[x_all.index(x)] for x in x_all if x in x_all_cellCycleCommitted]

#DZ
mcs_DZ_cellCycleCommitted = [int(a[idx_mcs]) for a in data_all_tuple if float(a[idx_cellCycleCommitted])==1 and float(a[idx_xCOM])<DZxCOM] #Extract all mcs points when cellCycleCommitted=1
mcs_DZ_cellCycleCommitted_regular = [i for i in mcs_DZ_cellCycleCommitted if i in mcs_regular_saves]
counts_DZ_cellCycleCommitted_at_mcs_regular = Counter(mcs_DZ_cellCycleCommitted_regular)
counts_DZ_cellCycleCommitted_at_mcs_regular = sorted(counts_DZ_cellCycleCommitted_at_mcs_regular.items()) #Sort by the key (mcs) and convert to a list of tuple
x_DZ_cellCycleCommitted = [tple[0] for tple in counts_DZ_cellCycleCommitted_at_mcs_regular]
y_DZ_cellCycleCommitted = [tple[1] for tple in counts_DZ_cellCycleCommitted_at_mcs_regular]
y_DZ_allcells = [y_DZ[x_DZ.index(x)] for x in x_DZ if x in x_DZ_cellCycleCommitted]

#LZ
mcs_LZ_cellCycleCommitted = [int(a[idx_mcs]) for a in data_all_tuple if float(a[idx_cellCycleCommitted])==1 and float(a[idx_xCOM])>LZxCOM] #Extract all mcs points when cellCycleCommitted=1
mcs_LZ_cellCycleCommitted_regular = [i for i in mcs_LZ_cellCycleCommitted if i in mcs_regular_saves]
counts_LZ_cellCycleCommitted_at_mcs_regular = Counter(mcs_LZ_cellCycleCommitted_regular)
counts_LZ_cellCycleCommitted_at_mcs_regular = sorted(counts_LZ_cellCycleCommitted_at_mcs_regular.items()) #Sort by the key (mcs) and convert to a list of tuple
x_LZ_cellCycleCommitted = [tple[0] for tple in counts_LZ_cellCycleCommitted_at_mcs_regular]
y_LZ_cellCycleCommitted = [tple[1] for tple in counts_LZ_cellCycleCommitted_at_mcs_regular]
y_LZ_allcells = [y_LZ[x_LZ.index(x)] for x in x_LZ if x in x_LZ_cellCycleCommitted]

# Counts
fig, ax = plt.subplots()
plt.plot(np.divide(x_all_cellCycleCommitted,100), np.divide(y_all_cellCycleCommitted,1000), color="blue", label="Total")
plt.plot(np.divide(x_DZ_cellCycleCommitted,100), np.divide(y_DZ_cellCycleCommitted,1000), color="red", label="DZ")
plt.plot(np.divide(x_LZ_cellCycleCommitted,100), np.divide(y_LZ_cellCycleCommitted,1000), color="green", label="LZ")
plt.xlabel("Time (x100 mcs)",fontsize=16)
plt.ylabel("Number of B cells",fontsize=16)
plt.title("B cells engaged in cell cycle",fontsize=16)
plt.legend(loc='upper left')
plt.xticks(range(0,721,240))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 3B.png', dpi=600)
plt.show()

# Percentage
fig, ax = plt.subplots()
plt.plot(np.divide(x_all_cellCycleCommitted,100), np.divide(y_all_cellCycleCommitted, y_all_allcells), color="blue", label="Total")
plt.plot(np.divide(x_DZ_cellCycleCommitted,100), np.divide(y_DZ_cellCycleCommitted, y_DZ_allcells), color="red", label="DZ")
plt.plot(np.divide(x_LZ_cellCycleCommitted,100), np.divide(y_LZ_cellCycleCommitted, y_LZ_allcells), color="green", label="LZ")
plt.xlabel("Time (x100 mcs)",fontsize=16)
plt.ylabel("Fraction of B cells",fontsize=16)
plt.title("B cells engaged in cell cycle",fontsize=16)
plt.legend(loc='upper right')
plt.xticks(range(0,721,240))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 3C.png', dpi=600)
plt.show()


"""
Fig. 3E, 3F - Cell birth
"""

# Obtain mcs at new births in DZ, LZ, and midZ
mcs_at_birth = []
mcs_at_birth_DZ = []
mcs_at_birth_midZ = []
mcs_at_birth_LZ = []
xCOM_at_birth = []

for file in cellIDs_alive_cytokinesis: #Each cytokinesis event creates 2 cells, but there is only one net cell number increase per cytokinesis event. So it is convenient to just use these cells for counting new born cells. 
    lines = readFile(file)
    lastLine = lines[-1]
    birthmcs = lastLine.split()[idx_mcs]
    birthxCOM = lastLine.split()[idx_xCOM]
    if is_float(birthmcs):
        mcs_at_birth.append(int(birthmcs))
        xCOM_at_birth.append(float(birthxCOM))
        if float(birthxCOM) < DZxCOM:
            mcs_at_birth_DZ.append(int(birthmcs))
        elif float(birthxCOM) > LZxCOM:
            mcs_at_birth_LZ.append(int(birthmcs))
        else:
            mcs_at_birth_midZ.append(int(birthmcs))

# Calculate all instantaneous and cumulative births
counts_birth_at_mcs = Counter(mcs_at_birth) #Obtain unique mcs and corresponding count
counts_birth_at_mcs = sorted(counts_birth_at_mcs.items()) #Sort by the key (mcs) and convert to a list of tuple
mcs_at_birth_unique = [tple[0] for tple in counts_birth_at_mcs]
counts_at_birth = [tple[1] for tple in counts_birth_at_mcs]
counts_at_birth_cumulative = np.cumsum(counts_at_birth)*2 # For cumulatitive, the net births need to be doubled since two daughter cells are born from each division.

# Plot number of all instantaneous birth events at a mcs
fig, ax = plt.subplots()
plt.plot(np.divide(mcs_at_birth_unique,100), counts_at_birth, c="black")
plt.xlabel("Time (x100 mcs)")
plt.ylabel("Number of B cells born")
plt.title("Number of born cells at a mcs time point")
plt.xticks(range(0,721,240))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.show()

# Calculate instantaneous and cumulative births in DZ
counts_birth_at_mcs_DZ = Counter(mcs_at_birth_DZ) #Obtain unique mcs and corresponding count
counts_birth_at_mcs_DZ = sorted(counts_birth_at_mcs_DZ.items()) #Sort by the key (mcs) and convert to a list of tuple
mcs_at_birth_DZ_unique = [tple[0] for tple in counts_birth_at_mcs_DZ]
counts_at_birth_DZ = [tple[1] for tple in counts_birth_at_mcs_DZ]
counts_at_birth_DZ_cumulative = np.cumsum(counts_at_birth_DZ)*2 # For cumulatitive, the net births need to be doubled since two daughter cells are born from each division.

# Calculate instantaneous and cumulative births in LZ
counts_birth_at_mcs_LZ = Counter(mcs_at_birth_LZ) #Obtain unique mcs and corresponding count
counts_birth_at_mcs_LZ = sorted(counts_birth_at_mcs_LZ.items()) #Sort by the key (mcs) and convert to a list of tuple
mcs_at_birth_LZ_unique = [tple[0] for tple in counts_birth_at_mcs_LZ]
counts_at_birth_LZ = [tple[1] for tple in counts_birth_at_mcs_LZ]
counts_at_birth_LZ_cumulative = np.cumsum(counts_at_birth_LZ)*2 # For cumulatitive, the net births need to be doubled since two daughter cells are born from each division.

# Calculate instantaneous and cumulative births in midZ
counts_birth_at_mcs_midZ = Counter(mcs_at_birth_midZ) #Obtain unique mcs and corresponding count
counts_birth_at_mcs_midZ = sorted(counts_birth_at_mcs_midZ.items()) #Sort by the key (mcs) and convert to a list of tuple
mcs_at_birth_midZ_unique = [tple[0] for tple in counts_birth_at_mcs_midZ]
counts_at_birth_midZ = [tple[1] for tple in counts_birth_at_mcs_midZ]
counts_at_birth_midZ_cumulative = np.cumsum(counts_at_birth_midZ)*2 # For cumulatitive, the net births need to be doubled since two daughter cells are born from each division.


# Plot number of cumulative births over time
fig, ax = plt.subplots()
plt.plot(np.divide(mcs_at_birth_unique,100), np.divide(counts_at_birth_cumulative,1000), c="Black", label="Birth - Total")
plt.plot(np.divide(mcs_at_birth_DZ_unique,100), np.divide(counts_at_birth_DZ_cumulative,1000), c="red", label="Birth - DZ")
plt.plot(np.divide(mcs_at_birth_LZ_unique,100), np.divide(counts_at_birth_LZ,1000), c="green", label="Birth - LZ")
#plt.plot(np.divide(mcs_at_birth_midZ_unique,100), np.divide(counts_at_birth_midZ,1000), c="gray", label="Birth - midZ")
plt.xlabel("Time (x100 mcs)", fontsize=16)
plt.ylabel("Cumulative births (x1000)", fontsize=16)
plt.title("Cumulative B cell birth", fontsize=16)
plt.legend(loc='upper left')
plt.xticks(range(0,721,240))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 3F.png', dpi=600)
plt.show()


# Plot births in unit time by using histogram
bin_size = save_interval * 40
fig, axL = plt.subplots()
axR = axL.twinx()  

countsAll, binsAll, bars = axL.hist(np.divide(mcs_at_birth,100), histtype='step', linewidth=1, edgecolor="black", bins=range(0,int(mcs_max/100),int(bin_size/100)), label="birth - Total")
countsDZ, binsDZ, bars = axL.hist(np.divide(mcs_at_birth_DZ,100), histtype='step', linewidth=1, edgecolor="red", bins=range(0,int(mcs_max/100),int(bin_size/100)), label="birth - DZ")
countsLZ, binsLZ, bars = axL.hist(np.divide(mcs_at_birth_LZ,100), histtype='step', linewidth=1, edgecolor="green", bins=range(0,int(mcs_max/100),int(bin_size/100)), label="birth - LZ")
#axL.hist(mcs_at_birth_midZ, histtype='step', linewidth=1, edgecolor="gray", bins=range(0,mcs_max,bin_size), label="birth - midZ")
axL.set_ylabel("Number of cell births per " + str(bin_size) + " mcs", fontsize=16)

# DZ_LZ_birth_ratio = np.divide(countsDZ, countsLZ)
# axR.stairs(np.divide(DZ_LZ_birth_ratio[1:],100), binsLZ[1:], ls=':', color="black", label="DZ:LZ birth ratio") #Skip the first one which ofteen produces na due to dividing by zero
# axR.set_yscale("log")
# axR.set_ylabel("DZ:LZ birth ratio")

hR, lR = axR.get_legend_handles_labels()
hL, lL = axL.get_legend_handles_labels()
axL.legend(hR+hL, lR+lL, loc='upper left', bbox_to_anchor=(0.0, 1.0))


axR.set_title("B cell birth rate", fontsize=16)
plt.xlim(0, binsAll[-1])
plt.xticks(range(0,720,240))
#ax.tick_params(axis='both', direction='in', labelsize=20)
axL.set_xlabel("Time (x100 mcs)", fontsize=16)
axL.tick_params(direction='in', labelsize=20)
axR.tick_params(direction='in', labelsize=20)
plt.savefig('Fig. 3E.png', dpi=600)
plt.show()







#%% Fig. 3H-3O
"""
Fig. 3I-3O - death counts
"""
# Obtain mcs at deaths for all deaths and DZ, LZ, midZ deaths
mcs_at_death = []
mcs_at_death_DZ = []
mcs_at_death_midZ = []
mcs_at_death_LZ = []
xCOM_at_death = []

for file in cellIDs_dead:
    lines = readFile(file)
    lastLine = lines[-1]
    deadmcs = lastLine.split()[idx_mcs]
    deadxCOM = lastLine.split()[idx_xCOM]
    if is_float(deadmcs):
        mcs_at_death.append(int(deadmcs))
        xCOM_at_death.append(float(deadxCOM))
        if float(deadxCOM) < DZxCOM:
            mcs_at_death_DZ.append(int(deadmcs))
        elif float(deadxCOM) > LZxCOM:
            mcs_at_death_LZ.append(int(deadmcs))
        else:
            mcs_at_death_midZ.append(int(deadmcs))

                
## All deaths   
counts_dead_at_mcs = Counter(mcs_at_death) #Obtain unique mcs and corresponding count
counts_dead_at_mcs = sorted(counts_dead_at_mcs.items()) #Sort by the key (mcs) and convert to a list of tuple
   
mcs_at_death_unique = [tple[0] for tple in counts_dead_at_mcs]
counts_at_death = [tple[1] for tple in counts_dead_at_mcs]
counts_at_death_cumulative = np.cumsum(counts_at_death)

# Calculate instantaneous death events by merging with zero count at most mcs at regular saves
mcs_regular_saves_with_zero_death = [i for i in mcs_regular_saves if i not in mcs_at_death_unique]
counts_regular_saves_with_zero_death = [0 for _ in range(len(mcs_regular_saves_with_zero_death))]

x_mcs = mcs_at_death_unique + mcs_regular_saves_with_zero_death
y_count_death = counts_at_death + counts_regular_saves_with_zero_death

zipped_list = zip(x_mcs, y_count_death)
zipped_sorted = sorted(zipped_list, key = lambda x: x[0]) #x[0] refers to x_mcs
x1_all, y1_all = zip(*zipped_sorted) #zip(*xxx) functions as unzip


## Cumulative death counts in DZ
counts_dead_at_mcs_DZ = Counter(mcs_at_death_DZ) #Obtain unique mcs and corresponding count
counts_dead_at_mcs_DZ = sorted(counts_dead_at_mcs_DZ.items()) #Sort by the key (mcs) and convert to a list of tuple
mcs_at_death_DZ_unique = [tple[0] for tple in counts_dead_at_mcs_DZ]
counts_at_death_DZ = [tple[1] for tple in counts_dead_at_mcs_DZ]
counts_at_death_DZ_cumulative = np.cumsum(counts_at_death_DZ)


## Cumulative death counts in LZ
counts_dead_at_mcs_LZ = Counter(mcs_at_death_LZ) #Obtain unique mcs and corresponding count
counts_dead_at_mcs_LZ = sorted(counts_dead_at_mcs_LZ.items()) #Sort by the key (mcs) and convert to a list of tuple
mcs_at_death_LZ_unique = [tple[0] for tple in counts_dead_at_mcs_LZ]
counts_at_death_LZ = [tple[1] for tple in counts_dead_at_mcs_LZ]
counts_at_death_LZ_cumulative = np.cumsum(counts_at_death_LZ)


## Cumulative death counts in midZ
counts_dead_at_mcs_midZ = Counter(mcs_at_death_midZ) #Obtain unique mcs and corresponding count
counts_dead_at_mcs_midZ = sorted(counts_dead_at_mcs_midZ.items()) #Sort by the key (mcs) and convert to a list of tuple
mcs_at_death_midZ_unique = [tple[0] for tple in counts_dead_at_mcs_midZ]
counts_at_death_midZ = [tple[1] for tple in counts_dead_at_mcs_midZ]
counts_at_death_midZ_cumulative = np.cumsum(counts_at_death_midZ)


# Plot number of instantaneous death events at a mcs
fig, ax = plt.subplots()
plt.plot(np.divide(x1_all,100), y1_all, c="blue")
plt.xlabel("Time (x100 mcs)", fontsize=16)
plt.ylabel("Number of B cell death", fontsize=16)
plt.title("Number of B cell death at a mcs", fontsize=16)
plt.legend(loc='upper left')
plt.xticks(range(0,721,240))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 3KS.png', dpi=600)
plt.show()





# Plot deaths in unit time by using histogram
bin_size = save_interval * 40
fig, axL = plt.subplots()
axR = axL.twinx()  

countsAll, binsAll, bars = axL.hist(np.divide(mcs_at_death,100), histtype='step', linewidth=1, edgecolor="blue", bins=range(0, int(mcs_max/100), int(bin_size/100)), label="Death - Total")
countsDZ, binsDZ, bars = axL.hist(np.divide(mcs_at_death_DZ,100), histtype='step', linewidth=1, edgecolor="red", bins=range(0, int(mcs_max/100), int(bin_size/100)), label="Death - DZ")
countsLZ, binsLZ, bars = axL.hist(np.divide(mcs_at_death_LZ,100), histtype='step', linewidth=1, edgecolor="green", bins=range(0, int(mcs_max/100), int(bin_size/100)), label="Death - LZ")
#axR.hist(mcs_at_death_midZ, histtype='step', linewidth=1, edgecolor="gray", bins=range(0,mcs_max,bin_size), label="Death - midZ")
axL.set_ylabel("B cell deaths/" + str(bin_size) + " mcs", fontsize=16)

#DZ_LZ_death_ratio = [DZ/LZ for DZ, LZ in zip(countsDZ, countsLZ)] #Old method to calculate DZ_LZ_death_ratio 
DZ_LZ_death_ratio = np.divide(countsDZ, countsLZ)
axR.stairs(DZ_LZ_death_ratio[1:], binsLZ[1:], ls=':', color="black", label="DZ:LZ death ratio") #Skip the first one which ofteen produces na due to dividing by zero
axR.set_yscale("log")
axR.set_ylabel("DZ:LZ death ratio", fontsize=16)

hR, lR = axR.get_legend_handles_labels()
hL, lL = axL.get_legend_handles_labels()
axR.legend(hR+hL, lR+lL, loc='upper left', bbox_to_anchor=(0.0, 1))

axL.set_xlabel("Time (x100 mcs)", fontsize=16)
axR.set_title("B cell death rate", fontsize=16)
plt.xlim(0, binsAll[-1])
plt.xticks(range(0,720,240))
axL.tick_params(direction='in', labelsize=20)
axR.tick_params(direction='in', labelsize=20)
plt.savefig('Fig. 3I.png', dpi=600)
plt.show()




# Plot fraction of cell deaths in unit time 

# Calculate the number of cells in a bin of mcs, note using the maximal cell number during each bin unit time as the population size for the denominator. Using mean or other single mcs-point values may lead to >100% death rate in bin_size mcs time.
countsAll_population = []
for a in binsAll[0:-1]*100:
    if a in x_all and a+bin_size in x_all:
        countsAll_population.append(max(y_all[x_all.index(a):x_all.index(a+bin_size)]))
    else:
        countsAll_population.append(0)
        
countsDZ_population = []
for a in binsDZ[0:-1]*100:
    if a in x_DZ and a+bin_size in x_DZ:
        countsDZ_population.append(max(y_DZ[x_DZ.index(a):x_DZ.index(a+bin_size)]))
    else:
        countsDZ_population.append(0)                                             
                                             
countsLZ_population = []
for a in binsLZ[0:-1]*100:
    if a in x_LZ and a+bin_size in x_LZ:
        countsLZ_population.append(max(y_LZ[x_LZ.index(a):x_LZ.index(a+bin_size)]))
    else:
        countsLZ_population.append(0)             

# Old method for calculating the number of cells in a bin of mcs, which only works when there are always cells in the bin of mcs in DZ or LZ.
# countsAll_population = [y_all[x_all.index(a)] for a in binsAll[0:-1] ] #Use the cell number at the beginning of each unit time as the population for the denominator
# countsDZ_population = [y_DZ[x_DZ.index(a)] for a in binsDZ[0:-1] ] #Use the cell number at the beginning of each unit time as the population for the denominator
# countsLZ_population = [y_LZ[x_LZ.index(a)] for a in binsLZ[0:-1] ] #Use the cell number at the beginning of each unit time as the population for the denominator
fig, ax = plt.subplots()
plt.stairs(np.divide(countsAll, countsAll_population), binsAll, color="blue", label="Death - Total")
plt.stairs(np.divide(countsDZ, countsDZ_population), binsDZ, color="red", label="Death - DZ")
plt.stairs(np.divide(countsLZ, countsLZ_population), binsLZ, color="green", label="Death - LZ")
plt.xlabel("Time (x100 mcs)", fontsize=16)
plt.ylabel("Fractional death/" + str(bin_size) + " mcs", fontsize=16)
plt.title("Fractional B cell death rate", fontsize=16)
plt.legend(loc='upper right')
plt.xlim(0, binsAll[-1])
plt.xticks(range(0,720,240))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 3J.png', dpi=600)
plt.show()


# Plot number of cumulative deaths over time
fig, ax = plt.subplots()
plt.plot(np.divide(mcs_at_death_unique,100), np.divide(counts_at_death_cumulative,1000), c="blue", label="Death - Total")
plt.plot(np.divide(mcs_at_death_DZ_unique,100), np.divide(counts_at_death_DZ_cumulative,1000), c="red", label="Death - DZ")
plt.plot(np.divide(mcs_at_death_LZ_unique,100), np.divide(counts_at_death_LZ_cumulative,1000), c="green", label="Death - LZ")
#plt.plot(mcs_at_death_midZ_unique, counts_at_death_midZ_cumulative, c="gray", label="Death - midZ")
plt.xlabel("Time (x100 mcs)", fontsize=16)
plt.ylabel("Cumulative deaths (x1000)", fontsize=16)
plt.title("Cumulative B cell death", fontsize=16)
plt.legend(loc='upper left')
plt.xticks(range(0,721,240))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 3K.png', dpi=600)
plt.show()




# Obtain mcs and xcom at deaths because of lethal mutation, negative selecltion, no access to FDC/T cells
# deathMutation
mcs_at_death_deathMutation = []
xCOM_at_death_deathMutation = []
for file in cellIDs_dead_deathMutation:
    lines = readFile(file)
    lastLine = lines[-1]
    deadmcs = lastLine.split()[idx_mcs]
    deadxCOM = lastLine.split()[idx_xCOM]
    if is_float(deadmcs):
        mcs_at_death_deathMutation.append(int(deadmcs))
        xCOM_at_death_deathMutation.append(float(deadxCOM))


# deathTimer
mcs_at_death_deathTimer = []
xCOM_at_death_deathTimer = []
for file in cellIDs_dead_deathTimer:
    lines = readFile(file)
    lastLine = lines[-1]
    deadmcs = lastLine.split()[idx_mcs]
    deadxCOM = lastLine.split()[idx_xCOM]
    if is_float(deadmcs):
        mcs_at_death_deathTimer.append(int(deadmcs))
        xCOM_at_death_deathTimer.append(float(deadxCOM))
        
        
# neg_selected
mcs_at_death_deathTimer_neg_selected = []
xCOM_at_death_deathTimer_neg_selected = []
for file in cellIDs_dead_deathTimer_neg_selected:
    lines = readFile(file)
    lastLine = lines[-1]
    deadmcs = lastLine.split()[idx_mcs]
    deadxCOM = lastLine.split()[idx_xCOM]
    if is_float(deadmcs):
        mcs_at_death_deathTimer_neg_selected.append(int(deadmcs))
        xCOM_at_death_deathTimer_neg_selected.append(float(deadxCOM))


# No access to T cells
mcs_at_death_deathTimer_noTtouch = []
xCOM_at_death_deathTimer_noTtouch = []
for file in cellIDs_dead_deathTimer_noTtouch:
    lines = readFile(file)
    lastLine = lines[-1]
    deadmcs = lastLine.split()[idx_mcs]
    deadxCOM = lastLine.split()[idx_xCOM]
    if is_float(deadmcs):
        mcs_at_death_deathTimer_noTtouch.append(int(deadmcs))
        xCOM_at_death_deathTimer_noTtouch.append(float(deadxCOM))


# Plot deaths in unit time by using histogram
bin_size = save_interval * 40
fig, ax = plt.subplots()
countsAll, bins, bars = plt.hist(np.divide(mcs_at_death,100), histtype='step', linewidth=1, edgecolor="black", bins=range(0, int(mcs_max/100), int(bin_size/100)), label="Death - Total")
countsDeathMutation, bins, bars = plt.hist(np.divide(mcs_at_death_deathMutation,100), histtype='step', linewidth=1, edgecolor="red", bins=range(0, int(mcs_max/100), int(bin_size/100)), label="Death - Lethal mutation")
countsDeathNegSelected, bins, bars = plt.hist(np.divide(mcs_at_death_deathTimer_neg_selected,100), histtype='step', linewidth=1, edgecolor="green", bins=range(0, int(mcs_max/100), int(bin_size/100)), label="Death - Neg selected")
countsDeathNoTtouch, bins, bars = plt.hist(np.divide(mcs_at_death_deathTimer_noTtouch,100), histtype='step', linewidth=1, edgecolor="blue", bins=range(0, int(mcs_max/100), int(bin_size/100)), label="Death - No T cell access")

plt.xlabel("Time (x100 mcs)", fontsize=16)
plt.ylabel("B cell deaths/" + str(bin_size) + " mcs", fontsize=16)
plt.title("B cell death rate", fontsize=16)
plt.xlim(0, binsAll[-1])
plt.xticks(range(0,720,240))
plt.legend( loc='upper left')
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 3L.png', dpi=600)
plt.show()




# Plot % death composition in unit time by using histogram
deathFractionStack = []
deathFractionStack.append(np.divide(countsDeathMutation, countsAll))
deathFractionStack.append(np.divide(countsDeathNegSelected, countsAll))
deathFractionStack.append(np.divide(countsDeathNoTtouch, countsAll))

fig, ax = plt.subplots()
plt.stackplot(bins[1:], deathFractionStack, colors = ["red", "green", "blue"], labels=["Death - Lethal mutation", "Death - Neg selected", "Death - No T cell access"])
plt.xlabel("Time (x100 mcs)", fontsize=16)
plt.ylabel("% of deaths in " + str(bin_size) + " mcs", fontsize=16)
plt.title("% composition of B cell death", fontsize=16)
plt.xlim(0, binsAll[-1])
plt.legend( loc='upper right')
plt.xticks(range(0,720,240))
plt.legend( loc='upper right')
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 3M.png', dpi=600)
plt.show()




# Scatter of xCOM and mcs when a cell dies for different reasons 
fig, ax = plt.subplots()
plt.scatter(xCOM_at_death_deathMutation, np.divide(mcs_at_death_deathMutation,100), s=4, label="Lethal mutation")
#plt.scatter(xCOM_at_death_deathTimer, mcs_at_death_deathTimer, s=8, label="Death timer expires")
plt.scatter(xCOM_at_death_deathTimer_noTtouch, np.divide(mcs_at_death_deathTimer_noTtouch,100), alpha=0.5, s=4, label="No access to T cells")
plt.scatter(xCOM_at_death_deathTimer_neg_selected, np.divide(mcs_at_death_deathTimer_neg_selected,100), alpha=0.3, s=4, label="Negative selected")
plt.xlabel("X-coordinate", fontsize=16)
plt.ylabel("Time (x100 mcs)", fontsize=16)
plt.title("Time and location at B cell death", fontsize=16)
plt.xticks(range(0,200 + xOffset*2+1,50))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.legend(loc='upper left')
plt.savefig('Fig. 3N.png', dpi=600)
plt.show()



# Histogram of xCOM when a cell dies for different reasons 
fig, ax = plt.subplots()
plt.hist(xCOM_at_death_deathMutation, histtype='stepfilled', alpha=1, edgecolor="black", bins=range(0,200 + xOffset*2,5),label="Lethal mutation")
plt.hist(xCOM_at_death_deathTimer_noTtouch, histtype='stepfilled', alpha=0.8, edgecolor="black", bins=range(0,200 + xOffset*2,5), label="No access to T cells")
plt.hist(xCOM_at_death_deathTimer_neg_selected, histtype='stepfilled', alpha=0.8, edgecolor="black", bins=range(0,200 + xOffset*2,5), label="Negative selected")
plt.xlabel('X-coordinate', fontsize=16)
plt.ylabel('Count', fontsize=16)
plt.title("Distribution of location at B cell death", fontsize=16)
plt.xticks(range(0,200 + xOffset*2+1,50))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.legend(loc='upper right')
plt.savefig('Fig. 3O.png', dpi=600)
plt.show()



# Scatter of xCOM and mcs when a cell dies for all reasons combined
plt.scatter(xCOM_at_death, mcs_at_death, s=8)
plt.xlabel("xCOM")
plt.ylabel("mcs")
plt.title("xCOM vs. mcs at death")
plt.xticks(range(0,200 + xOffset*2+1,50))
plt.show()


# Plot final volume of all cells and vanishing cells
volumes_last = []
volumes_last_DZ = []
volumes_last_LZ = []
xCOM_last = []
xCOM_last_DZ = []
xCOM_last_LZ = []
xCOM_at_vanish = []
volume_vanished = []
mcs_at_vanish = []
cellIDs_vanished = []

for file in fileTuple:
    lines = readFile(file)
    lastLine = lines[-1]
    lastLineList = lastLine.split()
    mcs = lastLineList[idx_mcs]
    xCOM = lastLineList[idx_xCOM]
    Alive = lastLineList[idx_Alive]
    Cytokinesis = lastLineList[idx_Cytokinesis]
    plasmaCell = lastLineList[idx_plasmaCell]
    actualVolume = lastLineList[idx_actualVolume]
    
    if is_float(mcs):
        volumes_last.append(float(actualVolume))
        xCOM_last.append(float(xCOM))
        
        if float(xCOM) < DZxCOM:
            volumes_last_DZ.append(float(actualVolume))
            xCOM_last_DZ.append(float(xCOM))
            
        if float(xCOM) > LZxCOM:
            volumes_last_LZ.append(float(actualVolume))
            xCOM_last_LZ.append(float(xCOM))
            
        
        if int(Alive) == 1 and int(Cytokinesis) == 0 and int(plasmaCell) == 0 and (int(mcs)+save_interval) <= mcs_max: # The last condition is to exclude cells whose last time point is the same as or right before the simulation stop time.
            xCOM_at_vanish.append(float(xCOM))
            volume_vanished.append(float(actualVolume))
            mcs_at_vanish.append(int(mcs))
            cellIDs_vanished.append(file)
            
         
# Scatter of final volume vs. xCOM of all cells
plt.scatter(xCOM_last, volumes_last, s=4)
plt.xlabel('xCOM')
plt.ylabel('Cell volume')
plt.title("Final volume vs. xCOM of all cells")
plt.show()

# Histogram of final volumes of all cells  
plt.hist(volumes_last, histtype='stepfilled', edgecolor="black", bins=range(0,int(max(volumes_last)+5),1))
plt.xlabel('Cell volume')
plt.ylabel('Frequency')
plt.title("Last volume of all cells")
plt.show()

"""
Fig. 3H - Histogram of final volumes of DZ vs LZ cells  
"""
fig, ax = plt.subplots()
plt.hist(volumes_last_DZ, histtype='stepfilled', edgecolor="black", bins=range(0,int(max(volumes_last)+5),1), label="DZ")
plt.hist(volumes_last_LZ, histtype='stepfilled', edgecolor="black",  alpha=0.5, bins=range(0,int(max(volumes_last)+5),1), label="LZ")
plt.xlabel('Cell volume', fontsize=16)
plt.ylabel('Count', fontsize=16)
plt.title("Last volume of B cells", fontsize=16)
plt.xticks(range(0,130,30))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.legend(loc='upper left')
plt.savefig('Fig. 3H.png', dpi=600)
plt.show()



# Scatter of mcs vs. xCOM at which cells vanished
plt.scatter(xCOM_at_vanish, mcs_at_vanish, s=5)
plt.xlabel('xCOM')
plt.ylabel('mcs')
plt.xticks(range(0, 200+xOffset*2+1, 50))
plt.title("xCOM vs. mcs of vanished cells")
plt.show()

# Histogram of final volume of cells that vanished
plt.hist(volume_vanished, bins=range(0,int(max(volume_vanished)+1),1))
plt.xlabel('Cell volume')
plt.ylabel('Frequency')
plt.title("Final volume of vanished cells")
plt.show()

# Scatter of final volume vs. xCOM of cells that vanished
plt.scatter(xCOM_at_vanish, volume_vanished, s=5)
plt.xlabel('xCOM')
plt.ylabel('Cell volume')
plt.title("Final volume vs. xCOM of vanished cells")
plt.show()

# Calculate and plot vanished cell counts
counts_vanish_at_mcs = Counter(mcs_at_vanish) #Obtain unique mcs and corresponding count
counts_vanish_at_mcs = sorted(counts_vanish_at_mcs.items()) #Sort by the key (mcs) and convert to a list of tuple
mcs_at_vanish_unique = [tple[0] for tple in counts_vanish_at_mcs]
counts_at_vanish = [tple[1] for tple in counts_vanish_at_mcs]
counts_at_vanish_cumulative = np.cumsum(counts_at_vanish)

bin_size = save_interval * 40
fig, axL = plt.subplots()
axR = axL.twinx()  

axR.plot(mcs_at_vanish_unique, counts_at_vanish_cumulative, label="Cumulative number of vanished cells")
axR.set_ylabel("Cumulative number of vanished cells")

countsvanish, bins, bars = axL.hist(mcs_at_vanish, histtype='step', linewidth=1, edgecolor="blue", bins=range(0,mcs_max,bin_size), label="Cells vanished / " + str(bin_size) + " mcs")
axL.set_ylabel("Cells vanished / " + str(bin_size) + " mcs")

hR, lR = axR.get_legend_handles_labels()
hL, lL = axL.get_legend_handles_labels()
axL.legend(hR+hL, lR+lL, loc='upper left')

axL.set_xlabel("Time (mcs)")
axL.set_title("Number of vanished cells")
plt.xlim(0, bins[-1])
plt.show()






#%% Fig. 3G, 4A-4C, S4B, S5B
"""
Fig. 4A-4C, S4B, S5B - Affinity and generation - Plot time-course affinity and generation of all cells at a regular save_interval, including average, standard deviation, 95% percentile, min and max.
"""

# Create a list containing mcs, affinityscores and Generation
mcs_affinity_generation_all = [[int(a[idx_mcs]), float(a[idx_affinityScore]), float(a[idx_Generation])] for a in data_all_tuple] #Extract all mcs points

# Extract affinity and generation of cells at regular save mcs
affinities_at_regular_saves = [[] for _ in range(len(mcs_regular_saves))] #Using * doesn't work because all sub-lists are just clones and are changed to same values simultaneously
generations_at_regular_saves = [[] for _ in range(len(mcs_regular_saves))] #Using * doesn't work because all sub-lists are just clones and are changed to same values simultaneously

for mcs, affinity, generation in mcs_affinity_generation_all:
    if mcs in mcs_regular_saves and affinity >= 0: #affinityScore = -1 is used to mark for the case of lethal mutation
        idx = mcs_regular_saves.index(mcs)
        affinities_at_regular_saves[idx].append(affinity)
    
    if mcs in mcs_regular_saves:
        idx = mcs_regular_saves.index(mcs)
        generations_at_regular_saves[idx].append(generation)

# Calculate statistics for affinity         
affinities_mean = [mean(i) for i in affinities_at_regular_saves]  
affinities_std = [np.std(i) for i in affinities_at_regular_saves]
affinities_p25 = [np.percentile(i,25) for i in affinities_at_regular_saves]
affinities_p75 = [np.percentile(i,75) for i in affinities_at_regular_saves]
affinities_p2dot5 = [np.percentile(i,2.5) for i in affinities_at_regular_saves]
affinities_p97dot5 = [np.percentile(i,97.5) for i in affinities_at_regular_saves]
affinities_min = [min(i) for i in affinities_at_regular_saves]  
affinities_max = [max(i) for i in affinities_at_regular_saves]   
affinities_initial = affinities_at_regular_saves[0]
affinities_final = affinities_at_regular_saves[-1]

# Calculate +/- std from mean of affinity    
affinities_plus_std = [xi + yi for xi, yi in zip(affinities_mean, affinities_std)]
affinities_minus_std = [xi - yi for xi, yi in zip(affinities_mean, affinities_std)]




# Plot affinity distributions over time
n = 10 #Plot every nth saved mcs time points
x = mcs_regular_saves
fig, ax = plt.subplots()
plt.plot(np.divide(x[::n],100), affinities_mean[::n], c="blue", label='Mean')
#plt.plot(x, affinities_plus_std[::n], c="red", ls = ':', label='Std')
#plt.plot(x, affinities_plus_std[::n], c="red", ls = ':')
plt.plot(np.divide(x[::n],100), affinities_p25[::n], c="red", ls = ':', label='25-75%')
plt.plot(np.divide(x[::n],100), affinities_p75[::n], c="red", ls = ':')
plt.plot(np.divide(x[::n],100), affinities_p2dot5[::n], c="green", ls = ':', label='2.5-97.5%')
plt.plot(np.divide(x[::n],100), affinities_p97dot5[::n], c="green", ls = ':')
plt.plot(np.divide(x[::n],100), affinities_min[::n], c="orange", ls = ':', label='min-max')
plt.plot(np.divide(x[::n],100), affinities_max[::n], c="orange", ls = ':')
plt.xlabel('Time (x100 mcs)', fontsize=16)
plt.ylabel('Affinity', fontsize=16)
plt.xticks(range(0,721,240))
plt.title("Affinity maturation", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.legend(loc='upper left')


# Plotting for plasma cells
mcs_at_plasmaCell = []
affinity_at_plasmaCell = []
for file in cellIDs_alive_plasmaCell:
    lines = readFile(file)
    lastLine = lines[-1]
    a = lastLine.split()[idx_mcs]
    b = lastLine.split()[idx_affinityScore]
    if is_float(a):
        mcs_at_plasmaCell.append(int(a))
        affinity_at_plasmaCell.append(float(b))

        
# Scatter plot of mcs and affinity when a cell differentiate into a plasma cell
plt.scatter(np.divide(mcs_at_plasmaCell,100), affinity_at_plasmaCell, s=4, c='gray', label='Plasma cells')
plt.xlabel("Time (x100 mcs)", fontsize=16)
plt.ylabel("Affinity", fontsize=16)
#plt.xlim(0, mcs_at_plasmaCell[-1])
#plt.ylim(3, max(affinity_at_plasmaCell)+1)
#plt.title("mcs vs. affinity at differentiation into plasma cells")
plt.xticks(range(0,721,240))
plt.yticks(np.arange(2.5,15,2.5))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.legend(loc='upper left')
plt.savefig('Fig. 4A.png', dpi=600) # Change file name to Fig. S4B or Fig. S5B as needed
plt.show()




# Plot overlay of plasma cell birth rate and cumulative number of plasam cells
counts_plasmaCell_at_mcs = Counter(mcs_at_plasmaCell) #Obtain unique mcs and corresponding count
counts_plasmaCell_at_mcs = sorted(counts_plasmaCell_at_mcs.items()) #Sort by the key (mcs) and convert to a list of tuple

mcs_at_plasmaCell_unique = [tple[0] for tple in counts_plasmaCell_at_mcs]
counts_at_plasmaCell = [tple[1] for tple in counts_plasmaCell_at_mcs]
counts_at_plasmaCell_cumulative = np.cumsum(counts_at_plasmaCell)

bin_size = save_interval * 40
fig, axL = plt.subplots()
axR = axL.twinx()  

axR.plot(np.divide(mcs_at_plasmaCell_unique,100), np.divide(counts_at_plasmaCell_cumulative,1000), label="Cumulative")
axR.set_ylabel("Cumulative PCs", fontsize=16)

countsPlasmaCell, bins, bars = axL.hist(np.divide(mcs_at_plasmaCell,100), histtype='step', linewidth=1, edgecolor="blue", bins=range(0, int(mcs_max/100), int(bin_size/100)), label="Production rate")
axL.set_ylabel("PCs produced/" + str(bin_size) + " mcs", fontsize=16)

hR, lR = axR.get_legend_handles_labels()
hL, lL = axL.get_legend_handles_labels()
axL.legend(hR+hL, lR+lL, loc='upper left')

axL.set_xlabel("Time (x100 mcs)", fontsize=16)
axL.set_title("Plasma cell production", fontsize=16)
plt.xlim(0, bins[-1])
plt.xticks(range(0,700,240))
axL.tick_params(direction='in', labelsize=20)
axR.tick_params(direction='in', labelsize=20)
plt.savefig('Fig. 4B.png', dpi=600)
plt.show()



# Plot fraction of plasma cells born in unit time out of entire population
fig, ax = plt.subplots()
countsAll_population = [y_all[x_all.index(a)] for a in binsAll[0:-1]*100 ] #Use the cell number at the beginning of each unit time as the population for the denominator
plt.stairs(np.divide(countsPlasmaCell, countsAll_population), binsAll, color="blue")
plt.xlabel("Time (x100 mcs)", fontsize=16)
plt.ylabel("Fractional PC production/" + str(bin_size) + " mcs", fontsize=16)
plt.title("Fractional PC production rate", fontsize=16)
plt.xlim(0, binsAll[-1])
plt.xticks(range(0,720,240))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 4C.png', dpi=600)
plt.show()



# Plot histogram of initial and final B cell antibody affinity
plt.hist(affinities_final, bins=np.arange(0,16,0.25), histtype='stepfilled', edgecolor="black", label='End of simulation')
plt.hist(affinities_initial, bins=np.arange(0,16,0.25), histtype='stepfilled', edgecolor="black", label='Initial')
plt.xlabel('Antibody Affinity')
plt.ylabel('Frequency')
plt.legend(loc='upper left')
plt.title('Antibody Affinity - Initial vs. End of simulation')
plt.show()



"""
Fig. 3G - cell generation
"""
# Calculate statistics for generations         
generations_mean = [mean(i) for i in generations_at_regular_saves]  
generations_std = [np.std(i) for i in generations_at_regular_saves]
generations_p25 = [np.percentile(i,25) for i in generations_at_regular_saves]
generations_p75 = [np.percentile(i,75) for i in generations_at_regular_saves]
generations_p2dot5 = [np.percentile(i,2.5) for i in generations_at_regular_saves]
generations_p97dot5 = [np.percentile(i,97.5) for i in generations_at_regular_saves]
generations_min = [min(i) for i in generations_at_regular_saves]  
generations_max = [max(i) for i in generations_at_regular_saves]   
generations_initial = generations_at_regular_saves[0]
generations_final = generations_at_regular_saves[-1]

# Calculate +/- std from mean of generation    
generations_plus_std = [xi + yi for xi, yi in zip(generations_mean, generations_std)]
generations_minus_std = [xi - yi for xi, yi in zip(generations_mean, generations_std)]

# Plot generation distributions over time
n = 10 #Plot every nth saved mcs time points
x = mcs_regular_saves
fig, ax = plt.subplots()
plt.plot(np.divide(x[::n],100), generations_mean[::n], c="blue", label='Mean')
#plt.plot(x, generations_plus_std[::n], c="red", ls = ':', label='Std')
#plt.plot(x, generations_plus_std[::n], c="red", ls = ':')
plt.plot(np.divide(x[::n],100), generations_p25[::n], c="red", ls = ':', label='25-75%')
plt.plot(np.divide(x[::n],100), generations_p75[::n], c="red", ls = ':')
plt.plot(np.divide(x[::n],100), generations_p2dot5[::n], c="green", ls = ':', label='2.5-97.5%')
plt.plot(np.divide(x[::n],100), generations_p97dot5[::n], c="green", ls = ':')
plt.plot(np.divide(x[::n],100), generations_min[::n], c="orange", ls = ':', label='min-max')
plt.plot(np.divide(x[::n],100), generations_max[::n], c="orange", ls = ':')
plt.xlabel('Time (x100 mcs)', fontsize=16)
plt.ylabel('B cell generation', fontsize=16)
plt.xticks(range(0,721,240))
plt.title('B cell generation', fontsize=16)
plt.legend(loc='upper left')
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 3G.png', dpi=600)
plt.show()





#%% Fig. 7, S2

"""
Generating Figure 7 - Molecular-level results: violin plot, in situ
"""

observeTime = 70000//save_interval*save_interval  # This gives multiples of save_interval in mcs
data_at_a_mcs = [row for row in data_all_tuple if int(row[idx_mcs]) == observeTime]

xCOM = [float(row[idx_xCOM]) for row in data_at_a_mcs]
yCOM = [float(row[idx_yCOM]) for row in data_at_a_mcs]

variables_to_display = ['pMHCII','AKT','FOXO1', 'CD40', 'cRel','MYC','AP4','CXCR4','CXCR5','RelA','Blimp1','caspase3', 'mTOR' ]
data = []
labels = []
color_face = []

# Automatically extract each variable's data into a list
for variable_name in variables_to_display:
    idx = header.index(variable_name)
    globals()[variable_name+'_DZ'] = [float(row[idx]) for row in data_at_a_mcs if float(row[idx_xCOM])<midxCOM]
    globals()[variable_name+'_LZ'] = [float(row[idx]) for row in data_at_a_mcs if float(row[idx_xCOM])>=midxCOM]
    globals()[variable_name] = [float(row[idx]) for row in data_at_a_mcs]
    
    data.append(globals()[variable_name+'_DZ'])
    data.append(globals()[variable_name+'_LZ'])

    labels.append(variable_name)
    labels.append('')
    
    color_face.append(colors[3])
    color_face.append(colors[2])

# Normalize AP4 and mTOR
AP4_max = max(AP4)
mTOR_max = max(mTOR)
idx = 2*variables_to_display.index('AP4')
data[idx] = list(np.divide(data[idx], AP4_max)*100)
data[idx+1] = list(np.divide(data[idx+1], AP4_max)*100)
idx = 2*variables_to_display.index('mTOR')
data[idx] = list(np.divide(data[idx], mTOR_max)*100)
data[idx+1] = list(np.divide(data[idx+1], mTOR_max)*100)

## Violin Plot of all variables together
fig, axs = plt.subplots()
sns.stripplot(data=data[0:22], jitter=True, color='black', alpha=1, s=1) # display data points
vplot = sns.violinplot(data=data[0:22], palette=color_face, inner=None, cut=0, bw_adjust=0.5)
plt.xticks(range(len(data[0:22])), labels[0:22])
axs.set_aspect(0.025)
plt.tick_params(direction='out', labelsize=7)
plt.savefig('Fig. 7A.png', dpi=600)
plt.show()

       

## Simulate in situ hybridization
# modified hsv in 256 color class
hsv_modified = cm.get_cmap('hsv', 256)
newcmp = ListedColormap(hsv_modified(np.linspace(0.7, 0, 256))) #create new hsv colormaps in range of 0 (red) to 0.7 (blue)
for variable_name in variables_to_display:
    fig, ax = plt.subplots()
    plt.scatter(xCOM, yCOM, c= globals()[variable_name], s= 6, cmap=newcmp)
    ax.set_aspect(0.8)
    plt.xlabel('xCOM')
    plt.ylabel('yCOM')
    plt.xlim([0,250])
    plt.ylim([0,200])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.yticks(range(0, 201, 50))
    cbar = plt.colorbar(location='top', shrink=0.823)
    for t in cbar.ax.get_yticklabels():
         t.set_fontsize(16)
    cbar.set_label(variable_name)
    plt.savefig('Fig. 7B-J-' + variable_name + '.png', dpi=600)



"""
Generating Figure S2 - scatter plot of pairs of molecules
"""
## Scatter plots for fast view
for i, variable1 in enumerate(variables_to_display):
    for variable2 in variables_to_display[i+1:]:
        
        x = globals()[variable1+'_LZ']
        y = globals()[variable2+'_LZ']
        plt.scatter(x, y, c=colors[2])
        
        x = globals()[variable1+'_DZ']
        y = globals()[variable2+'_DZ']
        plt.scatter(x, y, c=colors[3], alpha=0.8)
        
        plt.xlabel(variable1)
        plt.ylabel(variable2)        
        plt.show()

# Fig. S2 - Manually scatter plot to save. Change variable1 and variable 2 names for different pairs in FIg. S2
variable1 = 'AKT'
variable2 = 'FOXO1'

fig, ax = plt.subplots()
x = globals()[variable1+'_LZ']
y = globals()[variable2+'_LZ']
plt.scatter(x, y, c=colors[2])

x = globals()[variable1+'_DZ']
y = globals()[variable2+'_DZ']
plt.scatter(x, y, c=colors[3],alpha=0.8)


plt.xlabel(variable1, fontsize=16)
plt.ylabel(variable2, fontsize=16)
plt.title(variable2 + ' vs ' + variable1, fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. S2-' + variable1 + '_'  + variable2 + '.png', dpi=600)  
plt.show()




#%% Fig 5A
"""
Generating Fig 5A - Cell fate fraction in LZ
"""
# Percentage of cells that return from LZ to DZ, divide in LZ, die in LZ, differentiate into plasma cells, and remain in LZ.
DZx = 90 + xOffset
midx = 100 + xOffset
LZx = 110 + xOffset

cellIDs_entered_LZ = []
cellIDs_returned_DZ = []
cellIDs_divided_in_LZ = []
cellIDs_cyle_committed_not_divided_in_LZ = []
cellIDs_died_in_LZ = []
cellIDs_plasmaCell_in_LZ = []
cellIDs_remained_in_LZ = []
cellIDs_decision_made = []
cellIDs_pos_selected = []
cellIDs_neg_selected = []


for file in cellIDs_ready_to_move_to_LZ:
    
    LZEntryTime1 = 1e10 #Set these two entry times to ensure only if both exists and the former>latter will the cell be counted as entered LZ.
    LZEntryTime2 = 0
    
    DZReturnTime1 = 1e10
    DZReturnTime2 = 0 
        
    lines = readFile(file)
    data_cell = []
    if len(lines) > 0:
        for line in lines:
            data_cell.append(line.split())
        
        i = 0
        for lineList in data_cell:
            i += 1
            xCOM = float(lineList[idx_xCOM])
            mcs = int(lineList[idx_mcs])
            if xCOM >= DZx and xCOM < midx:
                LZEntryTime1 = mcs
                break
            
        for lineList in data_cell[i:]: # i: is used here to progressively skip the previous lines to speed up the calculation.
            i += 1
            xCOM = float(lineList[idx_xCOM])
            mcs = int(lineList[idx_mcs])
            if xCOM >= LZx:
                LZEntryTime2 = mcs
                break
            
        if LZEntryTime2 > LZEntryTime1:
            
            cellIDs_entered_LZ.append(file)
            
            cellCycleCommitted = float(data_cell[-1][idx_cellCycleCommitted])
            Cytokinesis = int(data_cell[-1][idx_Cytokinesis])
            Alive = int(data_cell[-1][idx_Alive])
            plasmaCell = int(data_cell[-1][idx_plasmaCell])
            selectionDecision = int(data_cell[-1][idx_selectionDecision])
            xCOM_last = float(data_cell[-1][idx_xCOM])
            
            
            if cellCycleCommitted == 1:

            # To calculate B cells that have committed to cell cycle and returned to DZ
            
                # Old method
                # for lineList in data_cell[i:]:
                #     i += 1
                #     xCOM = float(lineList[idx_xCOM])
                #     mcs = int(lineList[idx_mcs])
                #     if (mcs > LZEntryTime2) and (xCOM <= LZx) and (xCOM > midx):
                #         DZReturnTime1 = mcs
                #         break
                
                # for lineList in data_cell[i:]:
                #     i += 1
                #     xCOM = float(lineList[idx_xCOM])
                #     mcs = int(lineList[idx_mcs])
                #     if (mcs > LZEntryTime2) and (xCOM <= midx):
                #         DZReturnTime2 = mcs
                #         break
                
                # if DZReturnTime2 > DZReturnTime1: 
                #     cellIDs_returned_DZ.append(file)
                
                # New method
                if xCOM_last <= midx: 
                    cellIDs_returned_DZ.append(file)    
                    
            # To calculate B cells that have committed to cell cycle and divided in LZ
                if Cytokinesis == 1 and xCOM_last > midx: #Note midX is used to define that the cell is within LZ
                    cellIDs_divided_in_LZ.append(file)
                   
            # To calculate B cells that have committed to cell cycle, are not divided and still in LZ
                if Cytokinesis == 0 and xCOM_last > midx: #Note midX is used to define that the cell is within LZ
                    cellIDs_cyle_committed_not_divided_in_LZ.append(file)        
          
          
        # To calculate B cells that died in LZ
            if Alive == 0 and xCOM_last > midx: #Note midX is used to define that the cell is within LZ
               cellIDs_died_in_LZ.append(file)
               
        # To calculate B cells that differentiated into plasma cells
            if plasmaCell == 1 and xCOM_last > midx: #Note midX is used to define that the cell is within LZ
               cellIDs_plasmaCell_in_LZ.append(file)
        
        # To calculate B cells that didn't die in LZ, didn't divide in LZ, didn't differentiate into plasma cells, and didn't return to DZ
            if Alive == 1 and cellCycleCommitted == 0 and plasmaCell == 0 and xCOM_last > midx: #Note midX is used to define that the cell is within LZ
               cellIDs_remained_in_LZ.append(file)
               
        # To calculate B cells that a selection decision was made
            if selectionDecision == 1:
               cellIDs_decision_made.append(file)
               
        # To calculate positively selected B cells
            if cellCycleCommitted == 1:
               cellIDs_pos_selected.append(file)
               
        # To calculate negatively selected B cells
            if selectionDecision == 1 and cellCycleCommitted == 0:
               cellIDs_neg_selected.append(file)
            
Cell_count_balance = len(cellIDs_returned_DZ) + len(cellIDs_divided_in_LZ) + len(cellIDs_cyle_committed_not_divided_in_LZ) + len(cellIDs_died_in_LZ) + len(cellIDs_plasmaCell_in_LZ) + len(cellIDs_remained_in_LZ)        
if Cell_count_balance != len(cellIDs_entered_LZ):
    print("Cell count not balanced!")
    print("All cells add up to " + str(Cell_count_balance))
    print("Cells that entered LZ is " + str(len(cellIDs_entered_LZ)))
else:
    print("Cell count balanced!")


percentage_Death = len(cellIDs_died_in_LZ) / len(cellIDs_entered_LZ)
percentage_Return = len(cellIDs_returned_DZ) / len(cellIDs_entered_LZ)
percentage_Differentiation = len(cellIDs_plasmaCell_in_LZ) / len(cellIDs_entered_LZ)
percentage_Remain = (len(cellIDs_remained_in_LZ) + len(cellIDs_cyle_committed_not_divided_in_LZ) + len(cellIDs_divided_in_LZ)) / len(cellIDs_entered_LZ)
vars = ["Death", "DZ re-entry", "PC", "Remaining"]
percentages = [percentage_Death, percentage_Return, percentage_Differentiation, percentage_Remain]

fig, ax = plt.subplots()
plt.bar(vars, percentages, color=['orange', 'green', 'purple', 'blue'])
plt.title('B cell fates in LZ', fontsize=16)
plt.ylabel('Cell fate fraction', fontsize=16)
plt.yticks(np.arange(0,0.6,0.1))
ax.tick_params(axis='x', direction='in', labelsize=15)
ax.tick_params(axis='y', direction='in', labelsize=20)
plt.savefig('Fig. 5A.png', dpi=600)
plt.show()



#%% Fig 4D
"""
Generating Fig 4D - Affinity distribution
"""
# Obtain affinity of alive cells
affinities_alive = []
for file in cellIDs_alive:
    lines = readFile(file)
    lastLine = lines[-1]
    lastScore = lastLine.split()[idx_affinityScore]
    if is_float(lastScore) and float(lastScore) >=0:
        affinities_alive.append(float(lastScore))

# Obtain affinity of plasma cells
affinities_alive_plasmaCell = []
for file in cellIDs_alive_plasmaCell:
    lines = readFile(file)
    lastLine = lines[-1]
    lastScore = lastLine.split()[idx_affinityScore]
    if is_float(lastScore):
        affinities_alive_plasmaCell.append(float(lastScore))

# Obtain affinity of dead cells dying of deathTimer
affinities_dead_deathTimer = []
for file in cellIDs_dead_deathTimer:
    lines = readFile(file)
    lastLine = lines[-1]
    lastScore = lastLine.split()[idx_affinityScore]
    if is_float(lastScore):
        affinities_dead_deathTimer.append(float(lastScore))

# Obtain affinity of dead cells dying of deathTimer after negatively selected by touching FDC and T cells
affinities_dead_deathTimer_neg_selected = []
for file in cellIDs_dead_deathTimer_neg_selected:
    lines = readFile(file)
    lastLine = lines[-1]
    lastScore = lastLine.split()[idx_affinityScore]
    if is_float(lastScore):
        affinities_dead_deathTimer_neg_selected.append(float(lastScore))

# Obtain affinity of dead cells dying of deathTimer for no chance of touching FDC and T cells
affinities_dead_deathTimer_noTtouch = []
for file in cellIDs_dead_deathTimer_noTtouch:
    lines = readFile(file)
    lastLine = lines[-1]
    lastScore = lastLine.split()[idx_affinityScore]
    if is_float(lastScore):
        affinities_dead_deathTimer_noTtouch.append(float(lastScore))


# Obtain affinity of positively selected cells
affinities_alive_pos_selected = []
for file in cellIDs_pos_selected:
    lines = readFile(file)
    lastLine = lines[-1]
    lastScore = lastLine.split()[idx_affinityScore]
    if is_float(lastScore):
        affinities_alive_pos_selected.append(float(lastScore))

# Obtain affinity of dead cells dying of lethal mutation
affinities_dead_deathMutation = []
for file in cellIDs_dead_deathMutation:
    lines = readFile(file)
    lastLine = lines[-1]
    lastScore = lastLine.split()[idx_affinityScore]
    if is_float(lastScore):
        affinities_dead_deathMutation.append(float(lastScore))

# Plot histograms of affinities
#plt.hist(affinities_dead_deathTimer, bins=np.arange(0,16,0.25), label='Dead - death timer')
fig, ax = plt.subplots()
plt.hist(affinities_dead_deathTimer_noTtouch, bins=np.arange(0,16,0.25), histtype='stepfilled',  edgecolor="black", label='Dead - no access to T cells')
plt.hist(affinities_alive_pos_selected, bins=np.arange(0,16,0.25), histtype='stepfilled',  alpha=0.7, edgecolor="black", label='Alive - pos selected')
plt.hist(affinities_alive_plasmaCell, bins=np.arange(-1,16,0.25), histtype='stepfilled',  alpha=0.8, edgecolor="black", label='Alive - plasma cells')
plt.hist(affinities_dead_deathTimer_neg_selected, bins=np.arange(0,16,0.25), histtype='stepfilled',  alpha=0.8, edgecolor="black", label='Dead - neg selected')
#plt.hist(affinities_dead_deathMutation, bins=np.arange(-1,16,0.25), label='Dead - lethal mutation')
plt.xlabel('Affinity', fontsize=16)
plt.ylabel('Count', fontsize=16)
plt.xlim(0, 15)
#plt.xticks(range(0,720,240))
plt.title('Affinity distribution', fontsize=16)
plt.legend(loc='upper left')
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 4D.png', dpi=600)
plt.show()




#%% Fig. 5B-5O
"""
Generating Fig 5B-5O - Migration time and LZ residence time.
"""
# Calculate and plot DZ-to-LZ migration time
LZ_destination_xCOM = midxCOM
DZ_to_LZ_migrationStartTime = []
DZ_to_LZ_migrationEndTime = []
DZ_to_LZ_migrationTime = []
DZ_to_LZ_migrationStart_xCOM = []
cellIDs_DZtoLZ_migration = []

for file in cellIDs_entered_LZ:
    
    lines = readFile(file)
    data_cell = []
    if len(lines) > 0:
        for line in lines:
            data_cell.append(line.split())
        
        # To obtain the DZ_to_LZ_migrationStartTime
        firstLineList = data_cell[0]    # B cells starting to migrate to LZ will have cellCycleCommitted=0 on first row
        DZ_to_LZ_StartTime = -1
        if is_float(firstLineList[idx_mcs]):
            cellCycleCommitted = float(firstLineList[idx_cellCycleCommitted])
            mcs = int(firstLineList[idx_mcs])
            xCOM = float(firstLineList[idx_xCOM])
            if cellCycleCommitted == 0: # Not necessary, but good to have it checked again here.
                DZ_to_LZ_StartTime = mcs
                DZ_to_LZ_migrationStartTime.append(mcs)
                DZ_to_LZ_migrationStart_xCOM.append(xCOM)
                
        # To obtain the DZ_to_LZ_migrationEndTime 
        if DZ_to_LZ_StartTime > -1:
            for lineList in data_cell[1:]: # At least skip the first line
                if is_float(lineList[idx_mcs]):
                    xCOM = float(lineList[idx_xCOM])
                    mcs = int(lineList[idx_mcs])
                    if xCOM > LZ_destination_xCOM:
                        DZ_to_LZ_EndTime = mcs
                        DZ_to_LZ_migrationEndTime.append(mcs)
                        DZ_to_LZ_migrationTime.append(DZ_to_LZ_EndTime - DZ_to_LZ_StartTime)
                        cellIDs_DZtoLZ_migration.append(file)
                        break


mean_DZ_to_LZ_migrationTime = mean(DZ_to_LZ_migrationTime)
median_DZ_to_LZ_migrationTime = median(DZ_to_LZ_migrationTime)
std_DZ_to_LZ_migrationTime = stdev(DZ_to_LZ_migrationTime)

fig, ax = plt.subplots()
plt.hist(DZ_to_LZ_migrationTime, bins=range(0,int(max(DZ_to_LZ_migrationTime)+0),100), histtype='stepfilled',  alpha=1, color=colors[0], edgecolor="black")        
plt.xlabel('Migration time (mcs)', fontsize=16)
plt.ylabel('Count', fontsize=16)
plt.xticks(range(0,1300,400))
plt.title("DZ-to-LZ migration", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5B.png', dpi=600)
plt.show()

fig, ax = plt.subplots()
plt.scatter(DZ_to_LZ_migrationStart_xCOM, DZ_to_LZ_migrationTime, color=colors[0], s=4) 
#plt.tricontour(x, y, z, 15, linewidths=0.5, colors='k')
plt.xlabel('Migration start X-coordinate', fontsize=16)
plt.ylabel('Migration time (mcs)', fontsize=16)
plt.xlim(0, 250)
plt.yticks(range(0,1400,400))
plt.title("DZ-to-LZ migration", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5C.png', dpi=600)
plt.show()

fig, ax = plt.subplots()
plt.scatter(np.divide(np.subtract(DZ_to_LZ_migrationEndTime, DZ_to_LZ_migrationTime),100), DZ_to_LZ_migrationTime, color=colors[0], s=4) # Cannot use the start time since it may have more cells than the end time.
plt.xlabel('Migration start time (x100 mcs)', fontsize=16)
plt.ylabel('Migration time (mcs)', fontsize=16)
plt.xticks(range(0,721,240))
plt.yticks(range(0,1400,400))
plt.title("DZ-to-LZ migration", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5D.png', dpi=600)
plt.show()

fig, ax = plt.subplots()
plt.scatter(DZ_to_LZ_migrationStart_xCOM, np.divide(np.subtract(DZ_to_LZ_migrationEndTime, DZ_to_LZ_migrationTime),100), color=colors[0], s=4) 
plt.xlabel('Migration start X-coordinate', fontsize=16)
plt.ylabel('Migration start time (x100 mcs)', fontsize=16)
plt.xlim(0, 250)
plt.yticks(range(0,721,240))
plt.title("DZ-to-LZ migration", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5E.png', dpi=600)
plt.show()

fig, ax = plt.subplots()
plt.hist(DZ_to_LZ_migrationStart_xCOM, bins=range(0,int(max(DZ_to_LZ_migrationStart_xCOM)+20),1), histtype='stepfilled',  alpha=1, color=colors[0], edgecolor="black")        
plt.xlabel('DZ-to-LZ migration start location (x)', fontsize=16)
plt.ylabel('Frequency', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.title("DZ-to-LZ migration", fontsize=16)
plt.show()



# Calculate and plot LZ-to-DZ migration time
DZ_destination_xCOM = midxCOM
LZ_to_DZ_migrationStartTime = []
LZ_to_DZ_migrationEndTime = []
LZ_to_DZ_migrationTime = []
LZ_to_DZ_migrationStart_xCOM = []
cellIDs_LZtoDZ_migration = []
for file in cellIDs_returned_DZ:
    
    lines = readFile(file)
    data_cell = []
    if len(lines) > 0:
        for line in lines:
            data_cell.append(line.split())
        
        # To obtain the LZ_to_DZ_migrationStartTime 
        i = 0
        LZ_to_DZ_StartTime = -1
        for lineList in data_cell:
            if is_float(lineList[idx_mcs]):
                cellCycleCommitted = float(lineList[idx_cellCycleCommitted])
                mcs = int(lineList[idx_mcs])
                xCOM = float(lineList[idx_xCOM])
                if cellCycleCommitted == 1: # A new cell cycle starts when cellCycleCommitted first turns from 0 to 1 and CXCR4 starts to be induced.
                    LZ_to_DZ_StartTime = mcs
                    LZ_to_DZ_migrationStartTime.append(mcs)
                    LZ_to_DZ_migrationStart_xCOM.append(xCOM)
                    i = data_cell.index(lineList)
                    break
                
        # To obtain the LZ_to_DZ_migrationEndTime
        if LZ_to_DZ_StartTime > -1:
            for lineList in data_cell[i+1:]: # To start the scanning after LZ_to_DZ_migrationStartTime
                if is_float(lineList[idx_mcs]):
                    xCOM = float(lineList[idx_xCOM])
                    mcs = int(lineList[idx_mcs])
                    if xCOM < DZ_destination_xCOM:
                        LZ_to_DZ_EndTime = mcs
                        LZ_to_DZ_migrationEndTime.append(mcs)
                        LZ_to_DZ_migrationTime.append(LZ_to_DZ_EndTime - LZ_to_DZ_StartTime)
                        cellIDs_LZtoDZ_migration.append(file)
                        break


mean_LZ_to_DZ_migrationTime = mean(LZ_to_DZ_migrationTime)
median_LZ_to_DZ_migrationTime = median(LZ_to_DZ_migrationTime)
std_LZ_to_DZ_migrationTime = stdev(LZ_to_DZ_migrationTime)

fig, ax = plt.subplots()
plt.hist(LZ_to_DZ_migrationTime, bins=range(0,int(max(LZ_to_DZ_migrationTime)+0),100), histtype='stepfilled', alpha=1, color=colors[1], edgecolor="black")        
plt.xlabel('Migration time (mcs)', fontsize=16)
plt.ylabel('Cell count', fontsize=16)
plt.xticks(range(0,1300,400))
plt.title("LZ-to-DZ migration", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5F.png', dpi=600)
plt.show()

fig, ax = plt.subplots()
plt.scatter(LZ_to_DZ_migrationStart_xCOM, LZ_to_DZ_migrationTime, color=colors[1], s=4) 
plt.xlabel('Migration start X-coordinate', fontsize=16)
plt.ylabel('Migration time (mcs)', fontsize=16)
plt.xlim(0, 250)
plt.title("LZ-to-DZ migration", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5G.png', dpi=600)
plt.show()

fig, ax = plt.subplots()
plt.scatter(np.divide(np.subtract(LZ_to_DZ_migrationEndTime, LZ_to_DZ_migrationTime),100), LZ_to_DZ_migrationTime, color=colors[1], s=4) # Cannot use the start time since it may have more points than the end time.
plt.xlabel('Migration start time (x100 mcs)', fontsize=16)
plt.ylabel('Migration time (mcs)', fontsize=16)
plt.xticks(range(0,721,240))
plt.title("LZ-to-DZ migration", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5H.png', dpi=600)
plt.show()

fig, ax = plt.subplots()
plt.scatter(LZ_to_DZ_migrationStart_xCOM, np.divide(np.subtract(LZ_to_DZ_migrationEndTime, LZ_to_DZ_migrationTime),100), color=colors[1], s=4) 
plt.xlabel('Migration start X-coordinate', fontsize=16)
plt.ylabel('Migration start time (x100 mcs)', fontsize=16)
plt.xlim(0, 250)
plt.yticks(range(0,721,240))
plt.title("LZ_to_DZ migration", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5I.png', dpi=600)
plt.show()

fig, ax = plt.subplots()
plt.hist(LZ_to_DZ_migrationStart_xCOM, bins=range(0,int(max(LZ_to_DZ_migrationStart_xCOM)+20),1), histtype='stepfilled',  alpha=1, color=colors[1], edgecolor="black")        
plt.xlabel('LZ_to_DZ migration start location (x)', fontsize=16)
plt.ylabel('Frequency', fontsize=16)
plt.title("LZ_to_DZ migration", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.show()



# Calculate and plot LZ residence time for those B cells returning to DZ
xCOM_DZLZ_border = midxCOM
mcs_LZEntryTime = []
mcs_LZExitTime = []
LZResidenceTime = []
xCOM_LZFurthestReach = []
cellIDs_ResidenceTime = []
for file in cellIDs_returned_DZ:
    
    lines = readFile(file)
    data_cell = []
    if len(lines) > 0:
        for line in lines:
            data_cell.append(line.split())
        
        # To obtain the LZ entry time
        i = 0
        LZEntryTime = -1 #Reset to -1 to make sure there is an actual LZEntryTime for the current cell, not use the value from the previous cell.
        for lineList in data_cell:
            if is_float(lineList[idx_mcs]):
                xCOM = float(lineList[idx_xCOM])
                mcs = int(lineList[idx_mcs])
                if xCOM > xCOM_DZLZ_border:
                    LZEntryTime = mcs
                    mcs_LZEntryTime.append(mcs)
                    i = data_cell.index(lineList)
                    break
                
        # To obtain the LZ exit time
        if LZEntryTime > -1: #Use this to make sure there is an actual LZEntryTime for the current cell, not use the value from the previous cell.
            for lineList in data_cell[i+1:]: # To start the scanning after LZEntryTime
                if is_float(lineList[idx_mcs]):
                    xCOM = float(lineList[idx_xCOM])
                    mcs = int(lineList[idx_mcs])
                    if xCOM < xCOM_DZLZ_border and mcs > LZEntryTime + 100: #> at least greater than 100 mcs to avoid instantaneous pushed-back-force across the DZ-LZ border:
                        LZExitTime = mcs
                        mcs_LZExitTime.append(mcs)
                        LZResidenceTime.append(LZExitTime - LZEntryTime)
                        xCOM_LZFurthestReach.append(max([float(line[idx_xCOM]) for line in data_cell]))
                        cellIDs_ResidenceTime.append(file)
                        break
        
mean_LZResidenceTime = mean(LZResidenceTime)
median_LZResidenceTime = median(LZResidenceTime)
std_LZResidenceTime = stdev(LZResidenceTime)

fig, ax = plt.subplots()
plt.hist(LZResidenceTime, bins=range(0,int(max(LZResidenceTime)),100), histtype='stepfilled', color=colors[2], alpha=1, edgecolor="black")
plt.xlabel('LZ residence time (mcs)', fontsize=16)
plt.ylabel('Frequency', fontsize=16)
plt.xticks(range(0,2401,600))
plt.yticks(range(0,2401,600))
plt.title('LZ residence time', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5J.png', dpi=600)
plt.show()

fig, ax = plt.subplots()
plt.scatter(xCOM_LZFurthestReach, LZResidenceTime, color=colors[2], s=4) # Cannot use the entry time since it may have more points than the exit time.
plt.xlabel('LZ furthest reach X-coordinate', fontsize=16)
plt.ylabel('LZ residence time (mcs)', fontsize=16)
plt.xlim(0, 250)
plt.title(" LZ residence time vs LZ furthest reach", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5K.png', dpi=600)
plt.show()

fig, ax = plt.subplots()
plt.scatter(np.divide(np.subtract(mcs_LZExitTime, LZResidenceTime),100), LZResidenceTime, color=colors[2], s=4) # Cannot use the entry time since it may have more points than the exit time.
plt.xlabel('LZ entry time (mcs)', fontsize=16)
plt.ylabel('LZ residence time (mcs)', fontsize=16)
plt.title(" LZ residence time vs LZ entry time", fontsize=16)
plt.xticks(range(0,721,240))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5L.png', dpi=600)
plt.show()

fig, ax = plt.subplots()
plt.scatter(np.subtract(mcs_LZExitTime, LZResidenceTime), xCOM_LZFurthestReach, color=colors[2], s=4) # Cannot use the entry time since it may have more points than the exit time.
plt.xlabel('LZ entry time (mcs)', fontsize=16)
plt.ylabel('LZ furthest reach (x)', fontsize=16)
plt.title("LZ entry time vs LZ furthest reach", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.show()



# Calculate and plot DZ-to-LZ-to-DZ round trip time
DZ_destination_xCOM = midxCOM
mcs_RoundTripStart = []
mcs_RoundTripEnd = []
RoundTripTime = []
xCOM_RoundTripStart = []
xCOM_LZFurthestReach = []
cellIDs_RoundTrip = []
for file in cellIDs_returned_DZ:
    
    lines = readFile(file)
    data_cell = []
    if len(lines) > 0:
        for line in lines:
            data_cell.append(line.split())
        
        # To obtain the round trip start time and location
        firstLineList = data_cell[0]    # B cells starting to migrate to LZ will have cellCycleCommitted=0 on first row
        Round_Trip_StartTime = -1
        if is_float(firstLineList[idx_mcs]):
            cellCycleCommitted = float(firstLineList[idx_cellCycleCommitted])
            mcs = int(firstLineList[idx_mcs])
            xCOM = float(firstLineList[idx_xCOM])
            if cellCycleCommitted == 0: # Not necessary, but good to have it checked again here.
                Round_Trip_StartTime = mcs
                mcs_RoundTripStart.append(mcs)
                xCOM_RoundTripStart.append(xCOM)
                
        # To obtain the round trip end time
        if Round_Trip_StartTime > -1:
            for lineList in data_cell[1:]: # At least skip the first line
                if is_float(lineList[idx_mcs]):
                    cellCycleCommitted = float(lineList[idx_cellCycleCommitted])
                    xCOM = float(lineList[idx_xCOM])
                    mcs = int(lineList[idx_mcs])
                    if xCOM < DZ_destination_xCOM and mcs > Round_Trip_StartTime + 100 and cellCycleCommitted == 1: #> at least greater than 100 mcs to avoid instantaneous pushed-back-force across the DZ-LZ border. cellCycleCommitted ==  1 is necessary to ensure this is the returning cell after positively selected.
                        Round_Trip_EndTime = mcs
                        mcs_RoundTripEnd.append(mcs)
                        RoundTripTime.append(Round_Trip_EndTime - Round_Trip_StartTime)
                        xCOM_LZFurthestReach.append(max([float(line[idx_xCOM]) for line in data_cell]))
                        cellIDs_RoundTrip.append(file)
                        break

mean_RoundTripTime = mean(RoundTripTime)
median_RoundTripTime = median(RoundTripTime)
std_RoundTripTime = stdev(RoundTripTime)

fig, ax = plt.subplots()
plt.hist(RoundTripTime, bins=range(0,int(max(RoundTripTime)+200),100), histtype='stepfilled', color=colors[3], alpha=1, edgecolor="black")
plt.xlabel('Round trip time (mcs)', fontsize=16)
plt.ylabel('Cell count', fontsize=16)
plt.xticks(range(0,2401,800))
plt.title('DZ-LZ-DZ round trip time', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5M.png', dpi=600)
plt.show()

fig, ax = plt.subplots()
plt.scatter(np.subtract(RoundTripTime, LZResidenceTime), LZResidenceTime, color=colors[3], s=4)
plt.xlabel('DZ-to-LZ migration time (mcs)', fontsize=16)
plt.ylabel('Round Trip Time (mcs)', fontsize=16)
plt.title("DZ-to-LZ migration time vs DZ-LZ-DZ Round Trip Time", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.show()

fig, ax = plt.subplots()
plt.scatter(np.subtract(mcs_RoundTripEnd, RoundTripTime), RoundTripTime, color=colors[3], s=4) # Cannot use the entry time since it may have more points than the exit time.
plt.xlabel('Round trip start time (mcs)', fontsize=16)
plt.ylabel('DZ-LZ-DZ Round Trip Time (mcs)', fontsize=16)
plt.title("Round trip start time vs DZ-LZ-DZ Round Trip Time", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.show()

fig, ax = plt.subplots()
plt.scatter(xCOM_RoundTripStart, RoundTripTime, color=colors[3], s=4) # Cannot use the entry time since it may have more points than the exit time.
plt.xlabel('Round trip start location)', fontsize=16)
plt.ylabel('DZ-LZ-DZ Round Trip Time (mcs)', fontsize=16)
plt.title("Round trip start location vs DZ-LZ-DZ Round Trip Time", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.show()

fig, ax = plt.subplots()
plt.scatter(np.subtract(mcs_RoundTripEnd, RoundTripTime), xCOM_LZFurthestReach, color=colors[3], s=4) # Cannot use the entry time since it may have more points than the exit time.
plt.xlabel('Round trip start time (mcs)', fontsize=16)
plt.ylabel('LZ furthest reach', fontsize=16)
plt.title("Round trip start time vs LZ furthest reach", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.show()

fig, ax = plt.subplots()
plt.scatter(xCOM_LZFurthestReach, RoundTripTime, color=colors[3], s=4) # Cannot use the entry time since it may have more points than the exit time.
plt.xlabel('LZ furthest reach', fontsize=16)
plt.ylabel('Round trip Time (mcs)', fontsize=16)
plt.title("LZ furthest reach vs round trip time", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.show()


# Alternative simple method to calculate and plot DZ-to-LZ-to-DZ round trip time
# It should be the same compared with the one above and prefered
RoundTripTime_simple_method = [DZ_to_LZ_migrationTime[cellIDs_DZtoLZ_migration.index(cellID)] + LZResidenceTime[cellIDs_ResidenceTime.index(cellID)] for cellID in cellIDs_ResidenceTime if cellID in cellIDs_DZtoLZ_migration]
fig, ax = plt.subplots()
plt.hist(RoundTripTime_simple_method, bins=range(0,int(max(RoundTripTime_simple_method)+200),100), histtype='stepfilled', color=colors[3], alpha=1, edgecolor="black")
plt.xlabel('Round trip time (mcs)', fontsize=16)
plt.ylabel('Cell count', fontsize=16)
plt.xticks(range(0,2401,800))
plt.title('DZ-LZ-DZ round trip time - simple method', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5M_simple_method.png', dpi=600)
plt.show()

fig, ax = plt.subplots()
plt.scatter(RoundTripTime, RoundTripTime_simple_method, s=1)
plt.xlabel('Round Trip Time - simple method (mcs)', fontsize=16)
plt.ylabel('Round Trip Time (mcs)', fontsize=16)
plt.title("Round Trip Time time vs Round Trip Time -simple method", fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.show()



# Plot histograms of all transition time together
fig, ax = plt.subplots()
plt.hist(np.subtract(RoundTripTime_simple_method, LZResidenceTime), bins=range(0,int(max(LZResidenceTime)),100),histtype='stepfilled',  alpha=1, edgecolor="black", label = 'DZ-to-LZ migration time')
plt.hist(LZ_to_DZ_migrationTime, bins=range(0,int(max(LZResidenceTime)),100),histtype='stepfilled',  alpha=0.8, edgecolor="black", label = 'LZ-to-DZ migration time')
plt.hist(LZResidenceTime, bins=range(0,int(max(LZResidenceTime)),100), histtype='stepfilled',  alpha=0.7, edgecolor="black", label = 'LZ residence time')
plt.hist(RoundTripTime, bins=range(0,int(max(RoundTripTime)),100), histtype='stepfilled',  alpha=0.6, edgecolor="black", label = 'Round trip time')
plt.xlabel('Time (mcs)', fontsize=16)
plt.ylabel('Cell count', fontsize=16)
plt.xticks(range(0,2401,800))
plt.title('Transition time', fontsize=16)
plt.legend()
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5N.png', dpi=600)
plt.show()


# Calculat and plot the fractions of time spent in DZ-to-LZ migration, residence time, LZ-to-DZ migration
fraction_LZ_to_DZ_time = np.divide(LZ_to_DZ_migrationTime, RoundTripTime)
fraction_DZ_to_LZ_time = np.divide(np.subtract(RoundTripTime, LZResidenceTime), RoundTripTime) # Only need to subtract the residence time because it already includes the LZ-to-DZ migration time
fraction_LZ_residence_time = np.divide(LZResidenceTime, RoundTripTime)

mean_fraction_LZ_to_DZ_time = mean(fraction_LZ_to_DZ_time)
mean_fraction_DZ_to_LZ_time = mean(fraction_DZ_to_LZ_time)
mean_fraction_LZ_residence_time = mean(fraction_LZ_residence_time)


fig, ax = plt.subplots()
plt.hist(fraction_DZ_to_LZ_time, bins=np.arange(0,1,0.02), histtype='stepfilled',  alpha=1, edgecolor="black", label = 'DZ-to-LZ migration time')
plt.hist(fraction_LZ_to_DZ_time, bins=np.arange(0,1,0.02), histtype='stepfilled',  alpha=0.8, edgecolor="black", label = 'LZ-to-DZ migration time')
plt.hist(fraction_LZ_residence_time, bins=np.arange(0,1,0.02), histtype='stepfilled',  alpha=0.7, edgecolor="black", label = 'LZ residence time')
plt.xlabel('Fractionn of round trip time', fontsize=16)
plt.ylabel('Cell count', fontsize=16)
plt.title('Fractional composition of round trip time', fontsize=16)
plt.legend()
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 5O.png', dpi=600)
plt.show()




#%% Fig. 6A-6D, 3D
"""
Generating Fig 6A-6D - Calculate cell cycle length
"""
# Calculate and plot cell cycle lengths
cellCycleLengths_first = [] # First cell cycle in a proliferation burst
cellCycleLengths_subsequent = [] # Subsequent cell cycles in a proliferation burst
mcs_cellCycleStart_first = [] # Time when the first cell cycle in a proliferation burst is initiated
mcs_cellCycleStart_subsequent = [] # Time when where the subsequent cell cycles in a proliferation burst are initiated
xCOM_cellCycleStart_first = [] # Location where the first cell cycle in a proliferation burst is initiated
xCOM_cellCycleStart_subsequent = [] # Location where the subsequent cell cycles in a proliferation burst are initiated
mcs_cytokinesis_first = []
xCOM_cytokinesis_first = []
mcs_cytokinesis_subsequent = []
xCOM_cytokinesis_subsequent = []
for file in cellIDs_alive_cytokinesis:
    
    lines = readFile(file)
    data_cell = []
    if len(lines) > 0:
        for line in lines:
            data_cell.append(line.split())
        
        cellCycleEndTime = int(data_cell[-1][idx_mcs])
        cellCycleEndxCOM = float(data_cell[-1][idx_xCOM])
                
        for lineList in data_cell:
            if is_float(lineList[idx_mcs]):
                cellCycleCommitted = float(lineList[idx_cellCycleCommitted])
                mcs = int(lineList[idx_mcs])
                xCOM = float(lineList[idx_xCOM])
                if cellCycleCommitted == 1: # A new cell cycle starts when cellCycleCommitted first turns from 0 to 1.
                    cellCycleStartTime = mcs
                    if data_cell.index(lineList) == 0: #For cell cycles subsequent to the first one in a proliferation burst in which on the first line cellCycleCommitted is already 1.
                        cellCycleLengths_subsequent.append(cellCycleEndTime - cellCycleStartTime)
                        mcs_cellCycleStart_subsequent.append(mcs)
                        xCOM_cellCycleStart_subsequent.append(xCOM)
                        mcs_cytokinesis_subsequent.append(cellCycleEndTime)
                        xCOM_cytokinesis_subsequent.append(cellCycleEndxCOM)
                    else:
                        cellCycleLengths_first.append(cellCycleEndTime - cellCycleStartTime)
                        mcs_cellCycleStart_first.append(mcs)
                        xCOM_cellCycleStart_first.append(xCOM)
                        mcs_cytokinesis_first.append(cellCycleEndTime)
                        xCOM_cytokinesis_first.append(cellCycleEndxCOM)
                    break
        
cellCycleLengths = cellCycleLengths_first + cellCycleLengths_subsequent

mean_cellCycleLengths_subsequent = mean(cellCycleLengths_subsequent)
mean_cellCycleLengths_first = mean(cellCycleLengths_first)

std_cellCycleLengths_subsequent = stdev(cellCycleLengths_subsequent)
std_cellCycleLengths_first = stdev(cellCycleLengths_first)

# Plot xCOM at start vs. xCOM at cytokinesis
fig, ax = plt.subplots()
plt.scatter(xCOM_cytokinesis_first, xCOM_cellCycleStart_first,   s=2, label='First cycle')
plt.scatter(xCOM_cytokinesis_subsequent, xCOM_cellCycleStart_subsequent, s=2, label='Subsequent cycles')
plt.xlabel('Cytokenesis X-coordinate', fontsize=16)
plt.ylabel('Cell cycle start X-coordinate', fontsize=16)
plt.title("Locations at cell cycle start vs cytokinesis", fontsize=16)
plt.xticks(range(0, 200+xOffset*2+1, 50))
plt.yticks(range(0, 200+xOffset*2+1, 50))
plt.legend()
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 6A.png', dpi=600)
plt.show()

plt.hist(xCOM_cytokinesis_subsequent, bins=range(0,200 + xOffset*2 ,10), histtype='stepfilled', alpha=1, edgecolor="black", label='Subsequent cytokinesis')
plt.hist(xCOM_cytokinesis_first, bins=range(0,200 + xOffset*2 ,10), histtype='stepfilled', alpha=0.8, edgecolor="black", label='First cytokinesis')
plt.xlabel('xCOM', fontsize=16)
plt.ylabel('Frequency', fontsize=16)
plt.title("xCOM at cytokinesis in proliferation burst", fontsize=16)
plt.legend()
plt.show()

fig, ax = plt.subplots()
plt.hist(cellCycleLengths_first, bins=range(0,int(max(cellCycleLengths)+200),100), histtype='stepfilled', edgecolor="black", label='First cycle')   
plt.hist(cellCycleLengths_subsequent, bins=range(0,int(max(cellCycleLengths)+200),100), histtype='stepfilled', edgecolor="black", label='Subsequent cycles')        
plt.xlabel('Cell cycle length (mcs)', fontsize=16)
plt.ylabel('Cell count', fontsize=16)
plt.xticks(range(0, 2001, 500))
plt.yticks(range(0, 25001, 5000))
plt.title("Cell cycle length", fontsize=16)
plt.legend(loc='upper right')
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 6B.png', dpi=600)
plt.show()

plt.scatter(xCOM_cellCycleStart_subsequent, cellCycleLengths_subsequent, s=6, label='Subsequent cycle')
plt.scatter(xCOM_cellCycleStart_first, cellCycleLengths_first, s=6, label='First cycle')
plt.xlabel('Cell cycle start location (x)')
plt.ylabel('Cell cycle length (mcs)')
plt.title('Cell cycle length vs start location')
plt.legend(loc='upper left')
plt.show()

fig, ax = plt.subplots()
plt.scatter(np.divide(mcs_cellCycleStart_first,100), cellCycleLengths_first, s=6, label='First cycle')
plt.scatter(np.divide(mcs_cellCycleStart_subsequent,100), cellCycleLengths_subsequent, s=6, label='Subsequent cycle')
plt.xlabel('Cell cycle start time (mcs)', fontsize=16)
plt.ylabel('Cell cycle length (mcs)', fontsize=16)
plt.xticks(range(0, 721, 240))
plt.title('Cell cycle length vs start time', fontsize=16)
plt.legend(loc='upper left')
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 6C.png', dpi=600)
plt.show()

fig, ax = plt.subplots()
plt.scatter(np.divide(mcs_cellCycleStart_first,100), xCOM_cellCycleStart_first, s=6, label='First cycle')
plt.scatter(np.divide(mcs_cellCycleStart_subsequent,100), xCOM_cellCycleStart_subsequent, s=6, label='Subsequent cycle')
plt.xlabel('Cell cycle start time (mcs)', fontsize=16)
plt.ylabel('Cell cycle start X-coordinate', fontsize=16)
plt.xticks(range(0, 721, 240))
plt.yticks(range(0, 251, 50))
plt.title('Cell cycle start location vs start time', fontsize=16)
plt.legend(loc='best')
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 6D.png', dpi=600)
plt.show()


"""
Generating Fig 3D - mcs vs. xCOM at cytokinesis events
"""
fig, ax = plt.subplots()
plt.scatter(xCOM_cytokinesis_subsequent, np.divide(mcs_cytokinesis_subsequent,100), s=1, label='Subsequent cytokinesis')
plt.scatter(xCOM_cytokinesis_first, np.divide(mcs_cytokinesis_first,100),  s=1, alpha=0.8, label='First cytokinesis')
plt.xlabel('Cytokinesis location (x)', fontsize=16)
plt.ylabel('Cytokinesis time (x100 mcs)', fontsize=16)
plt.title("xCOM vs. mcs of cytokinesis")
plt.xticks(range(0, 200+xOffset*2+1, 50))
#plt.legend()
plt.yticks(range(0,721,240))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.legend(loc='upper right')
plt.savefig('Fig. 3D.png', dpi=600)
plt.show()



#%% Fig. 6E-6J, 8
"""
Fig. 6E-6H - Calculate burst size and its correlations with affinity
"""

proliferationBurstStarterCellIDs = cellIDs_returned_DZ # Use cell cycle committed in the middle in future to include those cells dividing in LZ

# Function to traverse binary tree till stopping at offersrping cells that is no longer cell cycle committed given a starter cell
def print_tree(start_id, id_list, ids):
    for id in id_list:
        if id.split('_')[1] == start_id and id.split('_')[1] == id.split('_')[2]:
            ids.append(id)
        if id.split('_')[1] == start_id and id.split('_')[1] != id.split('_')[2]:
            ids.append(id)
            
            # Extract cellCycleCommitted 
            data = []
            for line in readFile(id):
                data.append(line.split())
            x = [row[idx_cellCycleCommitted] for row in data]
            x = [num for num in x if is_float(num)]
            x = [float(num) for num in x]
            
            if x[0] == 1: # Continue to traverse only if the cell is still cell cycle committed
                print_tree(id.split('_')[2], id_list, ids)
    return ids


# Obtain cell list of proliferation bursts
burstList = []
n = 0
for file in proliferationBurstStarterCellIDs:
    print(n)
    starterID = file.split('_')[2]
    n += 1
    #print("************")
    cells_in_burst = print_tree(starterID, fileTuple, [])
    # if cells_in_burst == []: # If there are no offspring cells, skip to the next starterID
    #     continue
    if file not in cells_in_burst:
        cells_in_burst.insert(0, file)
    #print(cells_in_burst)
    
    burstList.append(cells_in_burst)
    

# Calculate average division numbers for each proliferation burst
averageDivisionNumbers = []
starterCellaffinityScores = []
starterCellMYCPeak = []
starterCellMYCAUC = []
starterCellAP4Last = []
starterCellAP4Peak = []
starterCellAP4AUC = []
starterCellIDs = []

for i in range(len(burstList)):

    # extract affinity scores of burst starter cells
    starterCell = readFile(burstList[i][0]) #Make sure in future [0] is always the starter cell
    starterCellIDs.append(burstList[i][0])
    if len(starterCell) > 0: # was > 1
        data = []
        for line in starterCell:
            data.append(line.split())
        
        x = [row[idx_affinityScore] for row in data]
        x = [num for num in x if is_float(num)]
        starterCellaffinityScores.append(float(x[-1]))
        
        y = [row[idx_MYC] for row in data]
        y = [num for num in y if is_float(num)]
        y = np.array(y, dtype=float)
        starterCellMYCPeak.append(max(y))
        starterCellMYCAUC.append(sum(y)*(len(data)-1)*save_interval)
        
        z = [row[idx_AP4] for row in data]
        z = [num for num in z if is_float(num)]
        z = np.array(z, dtype=float)
        starterCellAP4Last.append(z[-1]) # Last time point value of AP4
        starterCellAP4Peak.append(max(z))
        starterCellAP4AUC.append(sum(z)*(len(data)-1)*save_interval)
        


    cellsInABurst = burstList[i]
    #print(cellsInABurst)
    numCellsDivided = [0] * 150  # Number of cells that divided a certain number of times assuming maximum 150 divisions. For example, numCellsDivided[0] represents the number of cells that divided 0 times, and numCellsDivided[i] represents the number of cells that divided i times.

    motherIDs = []
    for fileName in cellsInABurst:
        mother = fileName.split('_')[1]
        motherIDs.append(mother)

    for fileName in cellsInABurst:
        if fileName.split('_')[2] not in motherIDs: # Record number and generation of terminal cells
            generation = fileName.split('_')[0]
            numCellsDivided[int(generation)] += 1
    
    averageDivisionNumber = 0
    generations = []
    for fileName in cellsInABurst:
        generations.append(int(fileName.split('_')[0]))
    starterCellGeneration = min(generations)

    for d in range(len(numCellsDivided)):
        averageDivisionNumber += (d-starterCellGeneration) * numCellsDivided[d] / (2 ** (d-starterCellGeneration))

    averageDivisionNumbers.append(averageDivisionNumber)


mean_averageDivisionNumbers = mean(averageDivisionNumbers)
median_averageDivisionNumbers = median(averageDivisionNumbers)
std_averageDivisionNumbers = stdev(averageDivisionNumbers)


# Group metrics into different affinity bins
affinity_vec = []
burstSize_vec = []
MYCPeak_vec = []
MYCAUC_vec = []
AP4Last_vec = []
AP4Peak_vec = []
AP4AUC_vec = []
LZResidenceTime_vec = []
for i in range(4, math.ceil(max(starterCellaffinityScores))+1):
    if i == 4:
        affinity = [x for x,y in zip(starterCellaffinityScores,averageDivisionNumbers) if  x < (i+1)]
        burstSize = [y for x,y in zip(starterCellaffinityScores,averageDivisionNumbers) if x < (i+1)]
        MYCPeak = [y for x,y in zip(starterCellaffinityScores,starterCellMYCPeak) if x < (i+1)]
        MYCAUC = [y for x,y in zip(starterCellaffinityScores,starterCellMYCAUC) if x < (i+1)]
        AP4Last = [y for x,y in zip(starterCellaffinityScores,starterCellAP4Last) if x < (i+1)]
        AP4Peak = [y for x,y in zip(starterCellaffinityScores,starterCellAP4Peak) if x < (i+1)]
        AP4AUC = [y for x,y in zip(starterCellaffinityScores,starterCellAP4AUC) if x < (i+1)]
        LZRT = [y for x,y in zip(starterCellaffinityScores,LZResidenceTime) if x < (i+1)]
        
    else:
        affinity = [x for x,y in zip(starterCellaffinityScores,averageDivisionNumbers) if x >=i and x < (i+1)]
        burstSize = [y for x,y in zip(starterCellaffinityScores,averageDivisionNumbers) if x >=i and x < (i+1)]
        MYCPeak = [y for x,y in zip(starterCellaffinityScores,starterCellMYCPeak) if x >=i and x < (i+1)]
        MYCAUC = [y for x,y in zip(starterCellaffinityScores,starterCellMYCAUC) if x >=i and x < (i+1)]
        AP4Last = [y for x,y in zip(starterCellaffinityScores,starterCellAP4Last) if x >=i and x < (i+1)]
        AP4Peak = [y for x,y in zip(starterCellaffinityScores,starterCellAP4Peak) if x >=i and x < (i+1)]
        AP4AUC = [y for x,y in zip(starterCellaffinityScores,starterCellAP4AUC) if x >=i and x < (i+1)]
        LZRT = [y for x,y in zip(starterCellaffinityScores,LZResidenceTime) if x >=i and x < (i+1)]
        
    #burstSize = [x for x in averageDivisionNumbers if starterCellaffinityScores[averageDivisionNumbers.index(x)] >=i and starterCellaffinityScores[averageDivisionNumbers.index(x)] < (i+1)+3]
    affinity_vec.append(affinity)
    burstSize_vec.append(burstSize)
    MYCPeak_vec.append(MYCPeak)
    MYCAUC_vec.append(MYCAUC)
    AP4Last_vec.append(AP4Last)
    AP4Peak_vec.append(AP4Peak)
    AP4AUC_vec.append(AP4AUC)
    LZResidenceTime_vec.append(LZRT)


# Calcculate for mean, std, sem for metrics in different affinity bins
mean_affinity = []
std_affinity = []
sem_affinity = []

top5perc_burstSize = []
mean_burstSize = []
std_burstSize = []
sem_burstSize = []

mean_MYCPeak = []
std_MYCPeak = []
sem_MYCPeak = []

mean_MYCAUC = []
std_MYCAUC = []
sem_MYCAUC = []

mean_AP4Peak = []
std_AP4Peak = []
sem_AP4Peak = []

mean_AP4AUC = []
std_AP4AUC = []
sem_AP4AUC = []

mean_AP4Last = []
std_AP4Last = []
sem_AP4Last = []

mean_LZRT = []
std_LZRT = []
sem_LZRT = []


for i in range(0,len(affinity_vec)):
    mean_affinity.append(mean(affinity_vec[i]))
    std_affinity.append(stdev(affinity_vec[i]))
    sem_affinity.append(stats.sem(affinity_vec[i]))
    
    top5perc_burstSize.append(np.percentile(burstSize_vec[i],95.0))
    mean_burstSize.append(mean(burstSize_vec[i]))
    std_burstSize.append(stdev(burstSize_vec[i]))
    sem_burstSize.append(stats.sem(burstSize_vec[i]))
    
    mean_MYCPeak.append(mean(MYCPeak_vec[i]))
    std_MYCPeak.append(stdev(MYCPeak_vec[i]))
    sem_MYCPeak.append(stats.sem(MYCPeak_vec[i]))
    
    mean_MYCAUC.append(mean(MYCAUC_vec[i]))
    std_MYCAUC.append(stdev(MYCAUC_vec[i]))
    sem_MYCAUC.append(stats.sem(MYCAUC_vec[i]))
    
    mean_AP4Peak.append(mean(AP4Peak_vec[i]))
    std_AP4Peak.append(stdev(AP4Peak_vec[i]))
    sem_AP4Peak.append(stats.sem(AP4Peak_vec[i]))
    
    mean_AP4AUC.append(mean(AP4AUC_vec[i]))
    std_AP4AUC.append(stdev(AP4AUC_vec[i]))
    sem_AP4AUC.append(stats.sem(AP4AUC_vec[i]))
    
    mean_AP4Last.append(mean(AP4Last_vec[i]))
    std_AP4Last.append(stdev(AP4Last_vec[i]))
    sem_AP4Last.append(stats.sem(AP4Last_vec[i]))
    
    mean_LZRT.append(mean(LZResidenceTime_vec[i]))
    std_LZRT.append(stdev(LZResidenceTime_vec[i]))
    sem_LZRT.append(stats.sem(LZResidenceTime_vec[i]))


## Burst size vs affinity
fig, ax = plt.subplots()
plt.hist(averageDivisionNumbers, edgecolor="black", bins=np.arange(0, math.ceil(max(averageDivisionNumbers)), 0.5)) #for floats numbers
#plt.hist(averageDivisionNumbers, bins=range(0, math.ceil(max(averageDivisionNumbers)+1), 0.5))
plt.xlabel('Burst size (division number)', fontsize=16)
plt.ylabel('Cell count', fontsize=16)
plt.title("Proliferative burst size", fontsize=16)
plt.xticks(range(0,6,1))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 6E.png', dpi=600)
plt.show()


fig, ax = plt.subplots()
plt.scatter(starterCellaffinityScores, averageDivisionNumbers, s=2)
plt.xlabel('Affinity', fontsize=16)
plt.ylabel('Burst size (division number)', fontsize=16)
plt.title("Proliferative burst size vs BCR affinity", fontsize=16)
plt.yticks(range(0,6,1))
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 6F.png', dpi=600)
plt.show()

# Violin plot
data = burstSize_vec
labels = np.floor(mean_affinity)
#sns.stripplot(data=data, jitter=True, color='black', alpha=0.2, s=3) # display data points
sns.violinplot(data=data, inner='box', cut=0) #bw_adjust=0.75 for smoothness of the plot
plt.xticks(range(len(data)), labels)
plt.xlabel('Affinity')
plt.ylabel('Burst size')
plt.title('Proliferative burst size vs BCR affinity')
plt.show()


# 95th % vs mean
fig, ax = plt.subplots()
plt.errorbar(mean_affinity, top5perc_burstSize, xerr=std_affinity, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('Affinity (mean +/- std)', fontsize=16)
plt.ylabel('Burst size at 95th percentile', fontsize=16)
plt.title('Proliferative burst size vs BCR affinity', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 6G.png', dpi=600)
plt.show()


# mean vs mean
fig, ax = plt.subplots()
plt.errorbar(mean_affinity, mean_burstSize, xerr=std_affinity, linestyle='None', marker='o', capsize=3)
plt.errorbar(mean_affinity, mean_burstSize, yerr=sem_burstSize, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('Affinity (mean +/- std)', fontsize=16)
plt.ylabel('Burst size (mean +/- sem)', fontsize=16)
plt.title('Proliferative burst size vs BCR affinity', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 6H.png', dpi=600)
plt.show()



"""
Fig. 8A-8D - Calculate relationships between MYC, AP4, affinity, and busrt size
"""
## MYC peak vs affinity
plt.scatter(starterCellaffinityScores, starterCellMYCPeak, s=2)
plt.xlabel('Affinity Scores')
plt.ylabel('MYC peak level')
plt.show()

data = MYCPeak_vec
labels = np.floor(mean_affinity)
#sns.stripplot(data=data, jitter=True, color='black', alpha=0.2, s=3) # display data points
fig, ax = plt.subplots()
sns.violinplot(data=data, cut=0)
plt.xticks(range(len(data)), labels)
plt.xlabel('Affinity', fontsize=16)
plt.ylabel('MYC peak level', fontsize=16)
plt.title('MYC peak level vs BCR affinity', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 8A.png', dpi=600)
plt.show()

plt.errorbar(mean_affinity, mean_MYCPeak, yerr=std_MYCPeak, linestyle='None', marker='o', capsize=3)
plt.errorbar(mean_affinity, mean_MYCPeak, xerr=std_affinity, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('Affinity (mean +/- std)')
plt.ylabel('MYC peak level (mean +/- std)')
plt.title('MYC peak level vs BCR affinity')
#plt.savefig('Fig. 6G.png', dpi=600)
plt.show()



## MYC AUC vs affinity
plt.scatter(starterCellaffinityScores, starterCellMYCAUC, s=2)
plt.xlabel('Affinity')
plt.ylabel('MYC AUC level')
plt.show()

data = MYCAUC_vec
labels = np.floor(mean_affinity)
fig, ax = plt.subplots()
#sns.stripplot(data=data, jitter=True, color='black', alpha=0.2, s=3) # display data points
sns.violinplot(data=data, cut=0)
plt.xticks(range(len(data)), labels)
plt.xlabel('Affinity', fontsize=16)
plt.ylabel('MYC AUC level', fontsize=16)
plt.title('MYC AUC level vs BCR affinity', fontsize=16)
ax.tick_params(axis='both',direction='in', labelsize=20)
plt.savefig('Fig. 8B.png', dpi=600)
plt.show()

plt.errorbar(mean_affinity, mean_MYCAUC, yerr=std_MYCAUC, linestyle='None', marker='o', capsize=3)
plt.errorbar(mean_affinity, mean_MYCAUC, xerr=std_affinity, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('Affinity (mean +/- std)')
plt.ylabel('MYC AUC level (mean +/- std)')
plt.title('Affinity vs. MYC AUC level')
plt.show()



# Burst size vs MYC peak level
plt.scatter(starterCellMYCPeak, averageDivisionNumbers, s=1)
plt.xlabel('MYC peak level')
plt.ylabel('Burst size')
plt.show()

plt.errorbar(mean_MYCPeak, mean_burstSize, yerr=sem_burstSize, linestyle='None', marker='o', capsize=3)
plt.errorbar(mean_MYCPeak, mean_burstSize, xerr=sem_MYCPeak, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('MYC peak level (mean +/- sem)')
plt.ylabel('Burst size (mean +/- sem)')
plt.title('MYC peak level vs. Burst size')
plt.show()

plt.errorbar(mean_MYCPeak, top5perc_burstSize, xerr=sem_MYCPeak, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('MYC peak level (mean +/- sem)')
plt.ylabel('Burst size at 95th percentile')
plt.title('MYC peak level vs. Burst size')
plt.show()



# Burst size vs MYC AUC level
plt.scatter(starterCellMYCAUC, averageDivisionNumbers, s=1)
plt.xlabel('MYC AUC level')
plt.ylabel('Burst size')
plt.show()

plt.errorbar(mean_MYCAUC, mean_burstSize, yerr=sem_burstSize, linestyle='None', marker='o', capsize=3)
plt.errorbar(mean_MYCAUC, mean_burstSize, xerr=sem_MYCAUC, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('MYC AUC level (mean +/- sem')
plt.ylabel('Burst size (mean +/- sem)')
plt.title('MYC AUC level vs. Burst size')
plt.show()

plt.errorbar(mean_MYCAUC, top5perc_burstSize, xerr=sem_MYCAUC, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('MYC AUC level (mean +/- sem)')
plt.ylabel('Burst size at 95th percentile')
plt.title('MYC AUC level vs. Burst size')
plt.show()



# Burst size vs AP4 peak level
plt.scatter(starterCellAP4Peak, averageDivisionNumbers, s=1)
plt.xlabel('AP4 peak level')
plt.ylabel('Burst size')
plt.show()

plt.errorbar(mean_AP4Peak, mean_burstSize, yerr=sem_burstSize, linestyle='None', marker='o', capsize=3)
plt.errorbar(mean_AP4Peak, mean_burstSize, xerr=sem_AP4Peak, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('AP4 peak level (mean +/- sem)')
plt.ylabel('Burst size (mean +/- sem)')
plt.title('AP4 peak level vs. Burst size')
plt.show()

fig, ax = plt.subplots()
plt.errorbar(mean_AP4Peak, top5perc_burstSize, xerr=sem_AP4Peak, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('AP4 peak level (mean +/- sem)', fontsize=16)
plt.ylabel('Burst size at 95th percentile', fontsize=16)
plt.title('Burst siz vs AP4 peak level', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 8D.png', dpi=600)
plt.show()



# Burst size vs AP4 AUC level
plt.scatter(starterCellAP4AUC, averageDivisionNumbers, s=1)
plt.xlabel('AP4 AUC level')
plt.ylabel('Burst size')
plt.show()

plt.errorbar(mean_AP4AUC, mean_burstSize, yerr=sem_burstSize, linestyle='None', marker='o', capsize=3)
plt.errorbar(mean_AP4AUC, mean_burstSize, xerr=sem_AP4AUC, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('AP4 AUC level (mean +/- sem)')
plt.ylabel('Burst size (mean +/- sem)')
plt.title('AP4 AUC level vs. Burst size')
plt.show()

plt.errorbar(mean_AP4AUC, top5perc_burstSize, xerr=sem_AP4AUC, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('AP4 AUC level (mean +/- sem)')
plt.ylabel('Burst size at 95th percentile')
plt.title('AP4 AUC level vs. Burst size')
plt.show()



# Burst size vs AP4 level at first cytokinesis
plt.scatter(starterCellAP4Last, averageDivisionNumbers, s=1)
plt.xlabel('AP4 level at first cytokinesis')
plt.ylabel('Burst size')
plt.show()

plt.errorbar(mean_AP4Last, mean_burstSize, yerr=sem_burstSize, linestyle='None', marker='o', capsize=3)
plt.errorbar(mean_AP4Last, mean_burstSize, xerr=sem_AP4Last, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('AP4 level at first cytokinesis (mean +/- sem)')
plt.ylabel('Burst size (mean +/- sem)')
plt.title('AP4 level at first cytokinesis vs. Burst size')
plt.show()

plt.errorbar(mean_AP4Last, top5perc_burstSize, xerr=sem_AP4Last, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('AP4 level at first cytokinesis (mean +/- sem)')
plt.ylabel('Burst size at 95th percentile')
plt.title('AP4 level at first cytokinesis vs. Burst size')
plt.show()

plt.scatter(starterCellMYCPeak, starterCellAP4Peak, s=1)
plt.xlabel('MYC Peak Levels')
plt.ylabel('AP4 Peak Levels')
plt.show()

fig, ax = plt.subplots()
plt.scatter(starterCellMYCAUC, starterCellAP4Peak, s=1)
plt.xlabel('MYC AUC level', fontsize=16)
plt.ylabel('AP4 peak level', fontsize=16)
plt.xticks(np.arange(0,6e6,1e6))
plt.title('AP4 vs MYC', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 8C.png', dpi=600)
plt.show()

plt.scatter(starterCellMYCAUC, starterCellAP4AUC, s=1)
plt.xlabel('MYC AUC Levels')
plt.ylabel('AP4 AUC levels')
plt.show()

plt.scatter(starterCellMYCAUC, starterCellAP4Last, s=1)
plt.xlabel('MYC AUC Levels')
plt.ylabel('AP4 Levels at first cytokinesis')
plt.show()



"""
Fig. 6I-6J - Calculate correlations between affinity and LZ residence time
"""
# LZ residence time vs affinity
fig, ax = plt.subplots()
plt.scatter(starterCellaffinityScores, LZResidenceTime, s=1)
plt.xlabel('Affinity', fontsize=16)
plt.ylabel('LZ residence time (mcs)', fontsize=16)
plt.title('LZ residence time vs BCR affinity', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 6I.png', dpi=600)
plt.show()

data = LZResidenceTime_vec
labels = np.floor(mean_affinity)
fig, ax = plt.subplots()
#sns.stripplot(data=data, jitter=True, color='black', alpha=0.2, s=3) # display data points
sns.violinplot(data=data, inner='box', cut=0) #bw_adjust=0.75 for smoothness of the plot
plt.xticks(range(len(data)), labels)
plt.xlabel('Affinity', fontsize=16)
plt.ylabel('LZ Residence Time (mcs)', fontsize=16)
plt.title('Affinity vs. LZ Residence Time', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.show()

fig, ax = plt.subplots()
plt.errorbar(mean_affinity, mean_LZRT, xerr=std_affinity, linestyle='None', marker='o', capsize=3)
plt.errorbar(mean_affinity, mean_LZRT, yerr=sem_LZRT, linestyle='None', markerfacecolor='black', markeredgecolor='black', marker='o', capsize=3)
plt.xlabel('Affinity (mean +/- std)', fontsize=16)
plt.ylabel('LZ residence time (mean +/- sem)', fontsize=16)
plt.title('LZ residence time vs BCR affinity', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 6J.png', dpi=600)
plt.show()



plt.scatter(LZResidenceTime, starterCellMYCPeak, s=1)
plt.xlabel('LZ Residence Time (mcs)')
plt.ylabel('MYC peak level')
plt.show()

plt.scatter(LZResidenceTime, starterCellMYCAUC, s=1)
plt.xlabel('LZ Residence Time (mcs)')
plt.ylabel('MYC AUC Level')
plt.show()

plt.scatter(LZResidenceTime, starterCellAP4Peak, s=1)
plt.xlabel('LZ Residence Time (mcs)')
plt.ylabel('AP4 peak level')
plt.show()

plt.scatter(LZResidenceTime, starterCellAP4AUC, s=1)
plt.xlabel('LZ Residence Time (mcs)')
plt.ylabel('AP4 AUC Level')
plt.show()

plt.scatter(LZResidenceTime, starterCellAP4Last, s=1)
plt.xlabel('LZ Residence Time (mcs)')
plt.ylabel('AP4 level at fist cytokinesis')
plt.show()


LZ_selection_time = np.subtract(LZResidenceTime, LZ_to_DZ_migrationTime)

plt.scatter(starterCellaffinityScores, LZ_selection_time, s=1)
plt.xlabel('Affinity')
plt.ylabel('LZ selection time (mcs)')
plt.show()

plt.scatter(LZ_selection_time, starterCellAP4Peak, s=1)
plt.xlabel('LZ Selection Time (mcs)')
plt.ylabel('AP4 peak level')
plt.show()

plt.scatter(LZ_selection_time, starterCellAP4AUC, s=1)
plt.xlabel('LZ Selection Time (mcs)')
plt.ylabel('AP4 AUC Level')
plt.show()

plt.scatter(LZ_selection_time, starterCellAP4Last, s=1)
plt.xlabel('LZ Selection Time (mcs)')
plt.ylabel('AP4 level at fist cytokinesis')
plt.show()


plt.scatter(np.subtract(mcs_LZExitTime, LZResidenceTime), LZ_selection_time, color=colors[2], s=4) # Cannot use the entry time since it may have more points than the exit time.
plt.xlabel('LZ entry time (mcs)')
plt.ylabel('LZ selection time (mcs)')
plt.title("LZ entry time vs LZ slection time")
#plt.savefig('Fig. 5L.png', dpi=600)
plt.show()




#%% Fig. 4E-4F, S1A-S1B, S1C-S1D, S4C-S4D, S5C-S5D
"""
Fig. 4E-4F, S1A-S1B, S1C-S1D, S4C-S4D, S5C-S5D - Calculating B cell clonality percentage composition
"""

# Extract seederID of cells as cloneID at regular save mcs
cloneIDs_at_regular_saves = [[] for _ in range(len(mcs_regular_saves))] #Beginning with an empty list. Using * doesn't work because all sub-lists are just clones and are changed to same values simultaneously
generations_at_regular_saves = [[] for _ in range(len(mcs_regular_saves))] #Beginning with an empty list. Using * doesn't work because all sub-lists are just clones and are changed to same values simultaneously

for line in data_all_tuple:
    mcs = int(line[idx_mcs])
    cloneID = int(line[idx_seederID])
    generation = int(line[idx_Generation])
    
    if mcs in mcs_regular_saves:
        idx = mcs_regular_saves.index(mcs)
        cloneIDs_at_regular_saves[idx].append(cloneID)
        generations_at_regular_saves[idx].append(generation)


# Calculate and stack plot clone size over time for each clone
cloneSizeStack = []
for seeder in cellIDs_seeder:
    cloneSize = []
    for i in range(len(mcs_regular_saves)):
        n = cloneIDs_at_regular_saves[i].count(int(seeder.split('_')[2])) # Count instances containing a specific seederID
        cloneSize.append(n)
    cloneSizeStack.append(cloneSize)

fig, ax = plt.subplots()
plt.stackplot(np.divide(mcs_regular_saves,100), cloneSizeStack)
plt.xlabel('Time (x100 mcs)', fontsize=16)
plt.ylabel('Clone size (B cell count)', fontsize=16)
plt.xticks(range(0,721,240))
plt.yticks(range(0,2001,500))
plt.title('B cell clone size', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 4E.png', dpi=600) #Change file name to Fig. S1A, Fig. S1C, Fig. S4C, Fig. S5C as needed
plt.show() 


# Calculate and stack plot clone fraction over time for each clone
cloneSizeTransposed = list(zip(*cloneSizeStack)) # Transpose from clone x mcs to mcs x clone
cloneFractionStack = []
for i in range(len(cloneSizeTransposed)):
    sumRow = sum(cloneSizeTransposed[i])

    fractionRow = []
    for j in cloneSizeTransposed[i]:
        if sumRow != 0:
            fractionRow.append(j/sumRow)
        else:
            fractionRow.append(0)
        
    cloneFractionStack.append(fractionRow)
cloneFractionStack = list(zip(*cloneFractionStack)) # Transpose from mcs x clone back to clone x mcs

fig, ax = plt.subplots()
plt.stackplot(np.divide(mcs_regular_saves,100), cloneFractionStack)
plt.xlabel('Time (x100 mcs)', fontsize=16)
plt.ylabel('Fraction of B cell clones', fontsize=16)
plt.xticks(range(0,721,240))
plt.title('Fractional composition of B cell clones', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.savefig('Fig. 4F.png', dpi=600) #Change file name to Fig. S1B, Fig. S1D, Fig. S4D, Fig. S5D as needed
plt.show()    


# Bar graph to show clone sizes at last mcs
cloneIDs_at_last_mcs = cloneIDs_at_regular_saves[-1] 
cloneIDs_at_last_mcs.sort()
cloneSizes_at_last_mcs = Counter(cloneIDs_at_last_mcs)
plt.bar(cloneSizes_at_last_mcs.keys(), cloneSizes_at_last_mcs.values())
# Old code: counts, bins, bars = plt.hist(cloneIDs_at_last_mcs, bins=range(min(cloneIDs_at_last_mcs), max(cloneIDs_at_last_mcs)+2, 1))
plt.xlabel('Clone ID', fontsize=16)
plt.ylabel('Number of offspring cells', fontsize=16)
ax.tick_params(axis='both', direction='in', labelsize=20)
plt.title("Clone size at last mcs", fontsize=16)
plt.show()  




#%% Plot time-course affinity of a selected clone at a regular save_interval, including average score, standard deviation, 95% percentile, min and max scores

seederID = "295" #31_cellData_V50 (190, 295, 277), #32_cellData_V50 (261, 239, 326)

# Function to traverse binary tree for a seeder cell
def clone_retrieve(start_id, id_list, ids):
    for id in id_list:
        if id.split('_')[1] == start_id and id.split('_')[1] == id.split('_')[2]:
            ids.append(id)
        
        if id.split('_')[1] == start_id and id.split('_')[1] != id.split('_')[2]:
            ids.append(id)
            
            clone_retrieve(id.split('_')[2], id_list, ids)
    return ids            
   
    
cellIDs_in_a_clone = clone_retrieve(seederID, fileTuple, [])

cellIDs_in_clone295 = cellIDs_in_a_clone

# Combine all data points for all cells in a clone
data_in_a_clone = []
for file in cellIDs_in_a_clone:
    for line in readFile(file):  # read in one line at a time from current file
        lineList = line.split()  # convert a line to a list of strings
        if is_float(lineList[idx_mcs]):  # check that the line is not a header (which is supposed to exist only in seeder cell files) by using mcs
            data_in_a_clone.append(lineList)

# Create a list containing mcs, affinityscores and Generation
mcs_affinity_generation_clone = [[int(a[idx_mcs]), float(a[idx_affinityScore]), float(a[idx_Generation])] for a in data_in_a_clone] #Extract all mcs points

# Extract affinity and generation of cells at regular save mcs
affinities_at_regular_saves = [[] for _ in range(len(mcs_regular_saves))] #Using * doesn't work because all sub-lists are just clones and are changed to same values simultaneously
generations_at_regular_saves = [[] for _ in range(len(mcs_regular_saves))] #Using * doesn't work because all sub-lists are just clones and are changed to same values simultaneously

for mcs, affinity, generation in mcs_affinity_generation_clone:
    if mcs in mcs_regular_saves and affinity >= 0: #affinityScore = -1 is used to mark for the case of lethal mutation
        idx = mcs_regular_saves.index(mcs)
        affinities_at_regular_saves[idx].append(affinity)
    
    if mcs in mcs_regular_saves:
        idx = mcs_regular_saves.index(mcs)
        generations_at_regular_saves[idx].append(generation)

# Calculate statistics for affinity         
affinities_mean = [mean(i) for i in affinities_at_regular_saves]  
affinities_std = [np.std(i) for i in affinities_at_regular_saves]
affinities_p25 = [np.percentile(i,25) for i in affinities_at_regular_saves]
affinities_p75 = [np.percentile(i,75) for i in affinities_at_regular_saves]
affinities_p2dot5 = [np.percentile(i,2.5) for i in affinities_at_regular_saves]
affinities_p97dot5 = [np.percentile(i,97.5) for i in affinities_at_regular_saves]
affinities_min = [min(i) for i in affinities_at_regular_saves]  
affinities_max = [max(i) for i in affinities_at_regular_saves]   
affinities_initial = affinities_at_regular_saves[0]
affinities_final = affinities_at_regular_saves[-1]

# Calculate +/- std from mean of affinity    
affinities_plus_std = [xi + yi for xi, yi in zip(affinities_mean, affinities_std)]
affinities_minus_std = [xi - yi for xi, yi in zip(affinities_mean, affinities_std)]

# Plot affinity distributions over time
n = 10 #Plot every nth saved mcs time points
x = mcs_regular_saves
plt.plot(x[::n], affinities_mean[::n], c="blue", label='Mean')
#plt.plot(x, affinities_plus_std[::n], c="red", ls = ':', label='Std')
#plt.plot(x, affinities_plus_std[::n], c="red", ls = ':')
plt.plot(x[::n], affinities_p25[::n], c="red", ls = ':', label='25-75%')
plt.plot(x[::n], affinities_p75[::n], c="red", ls = ':')
plt.plot(x[::n], affinities_p2dot5[::n], c="green", ls = ':', label='2.5-97.5%')
plt.plot(x[::n], affinities_p97dot5[::n], c="green", ls = ':')
plt.plot(x[::n], affinities_min[::n], c="orange", ls = ':', label='min-max')
plt.plot(x[::n], affinities_max[::n], c="orange", ls = ':')
plt.xlabel('Time (mcs)')
plt.ylabel('Antibody Affinity')
plt.legend(loc='upper left')
plt.title('Clone ID ' + seederID)
plt.show()


# Plotting for plasma cells
mcs_at_plasmaCell = []
affinity_at_plasmaCell = []
for file in cellIDs_alive_plasmaCell:
    lines = readFile(file)
    lastLine = lines[-1]
    a = lastLine.split()[idx_mcs]
    b = lastLine.split()[idx_affinityScore]
    if is_float(a):
        mcs_at_plasmaCell.append(int(a))
        affinity_at_plasmaCell.append(float(b))
   
# Scatter plot of mcs and affinity when a cell differentiate into a plasma cell
plt.scatter(mcs_at_plasmaCell, affinity_at_plasmaCell, s=8)
plt.xlabel("Time (mcs)")
plt.ylabel("Affinity")
plt.xlim(0, mcs_at_plasmaCell[-1])
plt.ylim(3, max(affinity_at_plasmaCell)+1)
plt.title('Clone ID ' + seederID + ': mcs vs. affinity at differentiation into plasma cells')
plt.show()

# Plot histogram of initial and final B cell antibody affinity
plt.hist(affinities_final, bins=np.arange(0,16,0.25), histtype='stepfilled', edgecolor="black", label='End of simulation')
#plt.hist(affinities_initial, bins=np.arange(0,16,0.25), histtype='stepfilled', edgecolor="black", label='Initial')
plt.xlabel('Antibody Affinity')
plt.ylabel('Frequency')
plt.legend(loc='upper left')
plt.title('Clone ID ' + seederID + ': Antibody Affinity - Initial vs. End of simulation')
plt.show()


# Plot histograms of final B cell antibody affinity of several top dominant clones
affinities_final_clone190 = affinities_final
affinities_final_clone295 = affinities_final
affinities_final_clone277 = affinities_final

plt.hist(affinities_final_clone190, bins=np.arange(0,16,0.25), histtype='stepfilled', color=colors[3], alpha=1, edgecolor="black", label='Clone ')
plt.hist(affinities_final_clone190, bins=np.arange(0,16,0.25), histtype='stepfilled', color=colors[3], alpha=1, edgecolor="black", label='Clone ')

plt.xlabel('Antibody Affinity')
plt.ylabel('Frequency')
plt.legend(loc='upper left')
plt.title('Antibody Affinity of Clones')
plt.show()



# To do:
#Calculate burst size distribution for a clone
#Calculate the positive selection fraction and DZ return fraction etc over time for fixed time intervals



# Check for duplicate cellID
cellIDss = []
for file in files:
    cellIDss.append(file.split('_')[2])
cellIDss_unique = set(cellIDss)

duplicates = [number for number in cellIDss if cellIDss.count(number) > 1]
unique_duplicates = list(set(duplicates))







#%% Additional figures
"""
Additional figures - Plot cell population turnover
"""
# Plot births, deaths, plasma cells, and vanished cells in unit time by using histogram
countsAllBirth, binsAllBirth, barsBirth = plt.hist(mcs_at_birth, histtype='step', linewidth=1, edgecolor="red", bins=range(0,mcs_max,bin_size), label="Birth") #The number could be overestimated since some daughter cells are not created due to vanishing effect.
countsAllDeath, binsAllDeath, barsDeath = plt.hist(mcs_at_death, histtype='step', linewidth=1, edgecolor="blue", bins=range(0,mcs_max,bin_size), label="Death")
countsAllPlasma, binsAllPlasma, barsAllPlasma = plt.hist(mcs_at_plasmaCell, histtype='step', linewidth=1, edgecolor="green", bins=range(0,mcs_max,bin_size), label="Plasma cells")
countsAllVanish, binsAllVanish, barsAllVanish = plt.hist(mcs_at_vanish, histtype='step', linewidth=1, edgecolor="orange", bins=range(0,mcs_max,bin_size), label="Cells vanished")
plt.xlabel("Time (mcs)")
plt.title("Number of cell turnover per " + str(bin_size) + " mcs")
plt.xlim(0, binsAllBirth[-1])
plt.legend()
plt.show()


# Plot cumulative births, deaths, plasma cells, and vanished cells
plt.plot(mcs_at_birth_unique, counts_at_birth_cumulative, color='red', label="Births") #The number could be overestimated due to use of doubling method since some daughter cells are not created due to vanishing effect.
plt.plot(mcs_at_death_unique, counts_at_death_cumulative, color='blue', label="Deaths")
plt.plot(mcs_at_plasmaCell_unique, counts_at_plasmaCell_cumulative, color='green', label="Plasma cells")
plt.plot(mcs_at_vanish_unique, counts_at_vanish_cumulative, color='orange', label="Cell vanished")
plt.xlabel("Time (mcs)")
plt.title("Cumulative number of cells")
plt.legend()
plt.show()




# Calculate cumulative DZ-to-LZ migration
counts_DZtoLZ_at_mcs = Counter(DZ_to_LZ_migrationEndTime) #Obtain unique mcs and corresponding count
counts_DZtoLZ_at_mcs = sorted(counts_DZtoLZ_at_mcs.items()) #Sort by the key (mcs) and convert to a list of tuple
mcs_at_DZtoLZ_unique = [tple[0] for tple in counts_DZtoLZ_at_mcs]
counts_at_DZtoLZ = [tple[1] for tple in counts_DZtoLZ_at_mcs]
counts_at_DZtoLZ_cumulative = np.cumsum(counts_at_DZtoLZ)

# Calculate cumulative LZ-to-DZ migration
counts_LZtoDZ_at_mcs = Counter(LZ_to_DZ_migrationEndTime) #Obtain unique mcs and corresponding count
counts_LZtoDZ_at_mcs = sorted(counts_LZtoDZ_at_mcs.items()) #Sort by the key (mcs) and convert to a list of tuple
mcs_at_LZtoDZ_unique = [tple[0] for tple in counts_LZtoDZ_at_mcs]
counts_at_LZtoDZ = [tple[1] for tple in counts_LZtoDZ_at_mcs]
counts_at_LZtoDZ_cumulative = np.cumsum(counts_at_LZtoDZ)


# Plot DZ turnover
countsAllBirthDZ, binsAllBirthDZ, barsBirthDZ = plt.hist(mcs_at_birth_DZ, histtype='step', linewidth=1, edgecolor="red", bins=range(0,mcs_max,bin_size), label="Birth")
countsAllDeathDZ, binsAllDeathDZ, barsDeathDZ = plt.hist(mcs_at_death_DZ, histtype='step', linewidth=1, edgecolor="blue", bins=range(0,mcs_max,bin_size), label="Death")
countsAllVanish, binsAllVanish, barsAllVanish = plt.hist(mcs_at_vanish, histtype='step', linewidth=1, edgecolor="orange", bins=range(0,mcs_max,bin_size), label="Cells vanished")
countsDZtoLZ, binsDZtoLZ, barsDZtoLZ = plt.hist(DZ_to_LZ_migrationEndTime, histtype='step', linewidth=1, edgecolor="green", bins=range(0,mcs_max,bin_size), label="DZ-to-LZ migration")
countsLZtoDZ, binsLZtoDZ, barsLZtoDZ = plt.hist(LZ_to_DZ_migrationEndTime, histtype='step', linewidth=1, edgecolor="black", bins=range(0,mcs_max,bin_size), label="LZ-to-DZ migration")
plt.xlabel("Time (mcs)")
plt.title("Number of DZ cell turnover per " + str(bin_size) + " mcs")
plt.xlim(0, binsAllBirth[-1])
plt.legend()
plt.savefig('Fig. 3G.png', dpi=600)
plt.show()

# Plot cumulative cells - this plot can be meaningless
plt.plot(mcs_at_birth_DZ_unique, counts_at_birth_DZ_cumulative, color='red', label="Births") #The number could be overestimated due to use of doubling method since some daughter cells are not created due to vanishing effect.
plt.plot(mcs_at_death_DZ_unique, counts_at_death_DZ_cumulative, color='blue', label="Deaths")
plt.plot(mcs_at_vanish_unique, counts_at_vanish_cumulative, color='orange', label="Cell vanished")
plt.plot(mcs_at_DZtoLZ_unique, counts_at_DZtoLZ_cumulative, color='green', label="DZ-to-LZ migration")
plt.plot(mcs_at_LZtoDZ_unique, counts_at_LZtoDZ_cumulative, color='black', label="LZ-to-DZ migration")
plt.xlabel("Time (mcs)")
plt.title("Cumulative number of DZ cells")
plt.legend()
plt.show()


# Plot LZ turnover
countsAllBirthLZ, binsAllBirthLZ, barsBirthLZ = plt.hist(mcs_at_birth_LZ, histtype='step', linewidth=1, edgecolor="red", bins=range(0,mcs_max,bin_size), label="Birth")
countsAllDeathLZ, binsAllDeathLZ, barsDeathLZ = plt.hist(mcs_at_death_LZ, histtype='step', linewidth=1, edgecolor="blue", bins=range(0,mcs_max,bin_size), label="Death")
countsAllPlasma, binsAllPlasma, barsAllPlasma = plt.hist(mcs_at_plasmaCell, histtype='step', linewidth=1, edgecolor="orange", bins=range(0,mcs_max,bin_size), label="Plasma cell")
countsDZtoLZ, binsDZtoLZ, barsDZtoLZ = plt.hist(DZ_to_LZ_migrationEndTime, histtype='step', linewidth=1, edgecolor="green", bins=range(0,mcs_max,bin_size), label="DZ-to-LZ migration")
countsLZtoDZ, binsLZtoDZ, barsLZtoDZ = plt.hist(LZ_to_DZ_migrationEndTime, histtype='step', linewidth=1, edgecolor="black", bins=range(0,mcs_max,bin_size), label="LZ-to-DZ migration")
plt.xlabel("Time (mcs)")
plt.title("Number of LZ cell turnover per " + str(bin_size) + " mcs")
plt.xlim(0, binsAllBirth[-1])
plt.legend()
plt.show()


# Plot cumulative cells - this plot can be meaningless
plt.plot(mcs_at_birth_LZ_unique, counts_at_birth_LZ_cumulative, color='red', label="Births") #The number could be overestimated due to use of doubling method since some daughter cells are not created due to vanishing effect.
plt.plot(mcs_at_death_LZ_unique, counts_at_death_LZ_cumulative, color='blue', label="Deaths")
plt.plot(mcs_at_plasmaCell_unique, counts_at_plasmaCell_cumulative, color='orange', label="Plasma cells")
plt.plot(mcs_at_DZtoLZ_unique, counts_at_DZtoLZ_cumulative, color='green', label="DZ-to-LZ migration")
plt.plot(mcs_at_LZtoDZ_unique, counts_at_LZtoDZ_cumulative, color='black', label="LZ-to-DZ migration")
plt.xlabel("Time (mcs)")
plt.title("Cumulative number of LZ cells")
plt.legend()
plt.show()


# Plot fraction returning to DZ, death, and plasma cell formation.
cellIDs_DZtoLZ_migration_unit_time = [[] for _ in range(len(binsDZtoLZ)-1)] #Beginning with an empty list. Using * doesn't work because all sub-lists are just clones and are changed to same values simultaneously
cellIDs_LZtoDZ_migration_unit_time = [[] for _ in range(len(binsDZtoLZ)-1)] 
cellIDs_Death_in_LZ_unit_time = [[] for _ in range(len(binsDZtoLZ)-1)]
cellIDs_dead_deathTimer_neg_selected_unit_time = [[] for _ in range(len(binsDZtoLZ)-1)]
cellIDs_dead_deathTimer_noTtouch_LZ_unit_time = [[] for _ in range(len(binsDZtoLZ)-1)]
cellIDs_Plasma_in_LZ_unit_time = [[] for _ in range(len(binsDZtoLZ)-1)]

binsDZtoLZ_list = list(binsDZtoLZ)
for t1 in DZ_to_LZ_migrationEndTime:
    for t2 in binsDZtoLZ_list[1:]:
        if t1 < t2:
            idx1 = binsDZtoLZ_list.index(t2)-1
            idx2 = DZ_to_LZ_migrationEndTime.index(t1)
            cellIDs_DZtoLZ_migration_unit_time[idx1].append(cellIDs_DZtoLZ_migration[idx2])
            break


fraction_DZ_return_unit_time =[]
fraction_Death_in_LZ_unit_time =[]
fraction_dead_deathTimer_neg_selected_unit_time = []
fraction_dead_deathTimer_noTtouch_LZ_unit_time = []
fraction_Plasma_in_LZ_unit_time = []
for i in range(len(cellIDs_DZtoLZ_migration_unit_time)):
    cellIDs_LZtoDZ_migration_unit_time[i] = [id for id in cellIDs_DZtoLZ_migration_unit_time[i] if id in cellIDs_LZtoDZ_migration]
    cellIDs_Death_in_LZ_unit_time[i] = [id for id in cellIDs_DZtoLZ_migration_unit_time[i] if id in cellIDs_died_in_LZ]
    cellIDs_dead_deathTimer_neg_selected_unit_time[i] = [id for id in cellIDs_DZtoLZ_migration_unit_time[i] if id in cellIDs_dead_deathTimer_neg_selected]
    cellIDs_dead_deathTimer_noTtouch_LZ_unit_time[i] = [id for id in cellIDs_DZtoLZ_migration_unit_time[i] if id in cellIDs_dead_deathTimer_noTtouch_LZ]
    cellIDs_Plasma_in_LZ_unit_time[i] = [id for id in cellIDs_DZtoLZ_migration_unit_time[i] if id in cellIDs_alive_plasmaCell]

    if len(cellIDs_DZtoLZ_migration_unit_time[i]) > 0:
        fraction_DZ_return_unit_time.append(len(cellIDs_LZtoDZ_migration_unit_time[i])/len(cellIDs_DZtoLZ_migration_unit_time[i]))
        fraction_Death_in_LZ_unit_time.append(len(cellIDs_Death_in_LZ_unit_time[i])/len(cellIDs_DZtoLZ_migration_unit_time[i]))
        fraction_dead_deathTimer_neg_selected_unit_time.append(len(cellIDs_dead_deathTimer_neg_selected_unit_time[i])/len(cellIDs_DZtoLZ_migration_unit_time[i]))
        fraction_dead_deathTimer_noTtouch_LZ_unit_time.append(len(cellIDs_dead_deathTimer_noTtouch_LZ_unit_time[i])/len(cellIDs_DZtoLZ_migration_unit_time[i]))
        fraction_Plasma_in_LZ_unit_time.append(len(cellIDs_Plasma_in_LZ_unit_time[i])/len(cellIDs_DZtoLZ_migration_unit_time[i]))

    else:
        fraction_DZ_return_unit_time.append(float('nan'))
        fraction_Death_in_LZ_unit_time.append(float('nan'))
        fraction_dead_deathTimer_neg_selected_unit_time.append(float('nan'))
        fraction_dead_deathTimer_noTtouch_LZ_unit_time.append(float('nan'))
        fraction_Plasma_in_LZ_unit_time.append(float('nan'))
        

plt.stairs(fraction_DZ_return_unit_time, binsDZtoLZ, color="red", label="DZ return") 
plt.stairs(fraction_Death_in_LZ_unit_time, binsDZtoLZ, color="blue", label="Death")
plt.stairs(fraction_dead_deathTimer_neg_selected_unit_time, binsDZtoLZ, ls=':', color="blue", label="Death - neg selected")
plt.stairs(fraction_dead_deathTimer_noTtouch_LZ_unit_time, binsDZtoLZ, ls=':', color="orange", label="Death - no T cell touch")
plt.stairs(fraction_Plasma_in_LZ_unit_time, binsDZtoLZ, color="green", label="Plasma cell") 
plt.xlabel("Time (mcs)")
plt.title("Fraction of LZ cell fates per " + str(bin_size) + " mcs")
plt.xlim(0, binsDZtoLZ[-2])
plt.legend()
plt.show()



