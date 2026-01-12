from io import StringIO
import random
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from scipy.stats import wasserstein_distance
from itertools import product
import csv
import os

numpy.set_printoptions(legacy='1.25')

FILE_NAME = "" # global variable to be filled by runFromFile()
# Find median size of blocks
def findMedian(total:int, blockCount):
    GENOME_SIZE = int(jumpgenomesize.get())

    if total == 0:
        return 0

    median = total / 2.0
    # If the total is odd
    if total % 2 == 1:
        for i in range(GENOME_SIZE):
            median -= blockCount[i]
            if median <= 0:
                return i + 1

    # Else the total is even
    for i in range(GENOME_SIZE):
        median -= blockCount[i]
        if median < 0.500000000001:
            i += 1
            if median > -0.499999999999:
                for j in range(GENOME_SIZE):
                    if blockCount[j] != 0:
                        return (j + 1 + i) / 2.0
            return i

    return -1

# Display results in GUI
def display_results(blockCount):
    genomeSize = len(blockCount)
    # Calculate the total amount of blocks
    totalBlocks = 0
    for i in range(genomeSize):
        totalBlocks += blockCount[i]
    #print(totalBlocks)

    # Find block frequency
    frequency = [0] * genomeSize
    if totalBlocks > 0:
        for i in range(genomeSize):
            frequency[i] = blockCount[i] / totalBlocks

    # Find average
    avg = 0.0
    if totalBlocks > 0:
        for i in range(genomeSize):
            avg += blockCount[i]*(i+1)
        avg /= totalBlocks

    for widget in frame_right.winfo_children():
        widget.destroy()

    # Labels for average and median
    label_avg = tk.Label(frame_right, text=f"Average Block Length: {avg:.2f}", font=("Arial", 12))
    label_avg.pack(pady=5)
    label_median = tk.Label(frame_right, text=f"Median Block Length: {findMedian(totalBlocks, blockCount)}", font=("Arial", 12))
    label_median.pack(pady=5)

    # Formula from preposition 4
    P = numpy.exp((-1) * float(edgelength.get()))
    Q = 1-P
    formulaDistribution = [0] * genomeSize
    Z = Q + 4*(P**2)*Q/(1+P)
    formulaDistribution[0] = (Q + 4*(P**2)*(Q**2))/Z
    for k in range(2, genomeSize):
        formulaDistribution[k-1] = (4 * (P**(2*k)) * Q**2)/Z
    formulaDistribution[genomeSize-1] = 1 - sum(formulaDistribution)
    # print("Formula distribution:")
    # print(formulaDistribution)

    label_wasserstein = tk.Label(frame_right, text=f"Wasserstein Distance: {compareDistributions(frequency, formulaDistribution)}", font=("Arial", 12))
    label_wasserstein.pack(pady=5)

    # Prepare data for histogram
    lengths = [(i+1) for i in range(genomeSize) if blockCount[i] > 0]
    valuesReal = [frequency[i] for i in range(genomeSize) if frequency[i] > 0]
    valuesSim = [formulaDistribution[i] for i in range(genomeSize) if frequency[i] > 0]

    # Plot with matplotlib
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(lengths, valuesSim, color='red', label='Simulated data', marker='s')   # Red line with squares
    ax.plot(lengths, valuesReal, color='blue', label='Induced data', marker='o')   # Blue line with circles
    ax.set_title("Synteny block distribution")
    ax.set_xlabel("Synteny Block Length")
    ax.set_ylabel("Frequency")

    # Embed plot into Tkinter
    canvas = FigureCanvasTkAgg(fig, master=frame_right)
    canvas.draw()
    canvas.get_tk_widget().pack()

# Print all blocks
def printBlocks(blocks):
    for block in blocks:
        print(str(block[0]) + ": " + str(block[1]) + "-" + str(block[1]+block[3]-1) + ", " + str(block[2]) + "-" + str(block[2]+block[3]-1))

# Jump model
def jump(genome, edgelength):
    GENOME_SIZE = int(jumpgenomesize.get())
    # Determine how many jumps should occur based on a poisson distribution
    jumpcount = numpy.random.poisson(edgelength*GENOME_SIZE)

    # Perform the jumps
    output = genome.copy()
    for i in range(jumpcount):
        # Determine who jumps and to where
        who = random.randint(0,GENOME_SIZE-1)
        where = random.randint(0,GENOME_SIZE-1)
        # Reroll if it jumped to the exact same spot
        while (output[where] == who):
            where = random.randint(0,GENOME_SIZE-1)
        output.insert(where,output.pop(who))

    #print("Total jumps: " + str(jumpcount) + " with edge length " + str(edgelength))
    return output

# Find synteny blocks between 2 genomes in the jump model
def findSyntenyJump(genome1, genome2):
    GENOME_SIZE = int(jumpgenomesize.get())
    blockCount = [0] * GENOME_SIZE
    # block = []
    # blocks = []

    # Iterate through the genomes
    for i in range(GENOME_SIZE):
        length = 0
        longestBlockLength = 0
        # block = []

        # Skip deleted genes
        while (i < GENOME_SIZE and genome1[i] == -1):
            i += 1
        if (i == GENOME_SIZE):
            break
        
        for j in range(GENOME_SIZE):
            # First gene in block found, find block length
            if i < GENOME_SIZE and genome1[i] != -1 and genome1[i] == genome2[j]:
                length = 1
                while (i+length < GENOME_SIZE and j+length < GENOME_SIZE and genome1[i+length] != -1 and genome1[i+length] == genome2[j+length]):
                    length += 1

                # If the current block is longer than the previously longest block, replace it
                if (length > longestBlockLength):
                    longestBlockLength = length
                    startPos = j
                j = j+length

        # Add the block to the list of blocks
        if (longestBlockLength > 0):
            for k in range(longestBlockLength):
                # block.append(genome1[i+k])
                genome1[i+k] = -1
                genome2[startPos+k] = -1
            # blocks.append([block, i, startPos, length]) 
            blockCount[length-1] += 1

    return blockCount

# Find synteny blocks between 2 real genomes from input
def findSyntenyReal(genome1, genome2):
    genome1size = len(genome1)
    genome2size = len(genome2)
    # block = []
    blocks = []

    # Because the genomes are circular, move the last gene in the first genome to position 0 until it is no longer the start of a block
    push = 0
    j = 0
    while j < genome2size:
        if push < min(genome1size,genome2size) and j > -1 and genome2[j] == genome1[-1] and genome2[(j+1)%genome2size] == genome1[0]:
            if push == 0:
                startPos2 = (j+1)%genome2size
            push += 1
            genome1.insert(0, genome1.pop())
            j -= 2
        j += 1

    # If the genomes are the same, return early
    if min(genome1size,genome2size) == push:
        blocks.append([push,0,genome1size-1,startPos2,(startPos2-1)%genome2size])
        return blocks

    # Iterate through the genomes to find synteny blocks
    blocksTemp = []
    addedBlock = True
    while addedBlock:
        addedBlock = False
        i = 0
        while i < genome1size:
            length = 0
            longestBlockLength = 0

            # Skip deleted genes
            while (i < genome1size and genome1[i] == -1):
                i += 1
            if (i == genome1size): # reached end of genome
                break
            
            # Compare with second genome (regular form)
            for j in range(genome2size):
                # First gene in block found, find block length
                if i < genome1size and genome1[i] != -1 and genome1[i] == genome2[j]:
                    length = 1
                    while (genome1[(i+length)%genome1size] != -1 and genome1[(i+length)%genome1size] == genome2[(j+length)%genome2size]):
                        length += 1

                    # If the current block is longer than the previous longest block, replace it
                    if (length > longestBlockLength):
                        inversed = False
                        longestBlockLength = length
                        startPos1 = i
                        startPos2 = j
                    j = j+length
            
            # Compare with second genome (inversed form)
            for j in range(genome2size-1, 0, -1):
                # First gene in block found, find block length
                if i < genome1size and genome1[i] != -1 and genome1[i] == genome2[j]:
                    length = 1
                    while (length < genome1size and genome1[(i+length)%genome1size] != -1 and genome1[(i+length)%genome1size] == genome2[(j-length)%genome2size]):
                        length += 1

                    # If the current block is longer than the previous longest block, replace it
                    if (length > longestBlockLength):
                        inversed = True
                        longestBlockLength = length
                        startPos1 = i
                        startPos2 = j
                    j = j-length

            # Add the longest block found to the list of blocks
            if (longestBlockLength > 0):
                addedBlock = True
                if (inversed):
                    blocksTemp.append([(startPos1-push)%genome1size, (startPos1+longestBlockLength-push-1),
                                        startPos2, (startPos2-longestBlockLength+1), -longestBlockLength])
                else:
                    blocksTemp.append([(startPos1-push)%genome1size, (startPos1+longestBlockLength-push-1),
                                        startPos2, (startPos2+longestBlockLength-1), longestBlockLength])
                i += longestBlockLength
            # If no block found, delete the gene (because no block will ever be found for it)
            else:
                genome1[i] = -1
                i += 1

        if addedBlock:
            blocksTemp = sorted(blocksTemp, key=lambda x: abs(x[0]), reverse=True)
            clearOverlaps(blocksTemp)
            deleteGenes(genome1, genome2, blocksTemp)
            # Move all the new blocks to the list of all blocks
            while len(blocksTemp) != 0:
                blocks.append(blocksTemp.pop(0))
        
    blocks = sorted(blocks, key=lambda x: abs(x[4]), reverse=True)

    # Increase all indexes by 1 to start at 1 instead of 0
    for i in range(len(blocks)):
        blocks[i][0] += 1
        blocks[i][1] += 1
        blocks[i][2] += 1
        blocks[i][3] += 1
        blocks[i][4] = abs(blocks[i][4])
    
    return blocks

# Remove any overlapping blocks
def clearOverlaps(blocks):
    i = 0
    while i < len(blocks):
        j = i+1
        while j < len(blocks):
            # Block j begins after block i and ends before block i or it ends after block i begins and ends before i ends
            if doesOverlap(blocks[i][0],blocks[i][1],blocks[j][0],blocks[j][1]) or \
               doesOverlap(blocks[i][2],blocks[i][3],blocks[j][2],blocks[j][3]):
                del blocks[j]
            else:
                j += 1
        i += 1

# Check if there is overlap between blocks i and j
def doesOverlap(starti, endi, startj, endj):
    return (starti <= max(startj,endj)) and (endi >= min(startj,endj))

    # # If the block is inversed
    # if jIsInversed:
    #     temp = startj
    #     startj = endj
    #     endj = temp

    # if iIsInversed:
    #     temp = starti
    #     starti = endi
    #     endi = temp 

    # # Both blocks wrap around (therefore both must overlap at 0)
    # if starti > endi and startj > endj:
    #     return True

    # # Circular case where only block i wraps around
    # if starti > endi:
    #     return doesOverlap(starti,genomesize-1,startj,endj,genomesize) or doesOverlap(0,endi,startj,endj,genomesize)
    
    # # Circular case where only block j wraps around
    # if startj > endj:
    #     return doesOverlap(starti,endi,startj,genomesize-1,genomesize) or doesOverlap(starti,endi,0,endj,genomesize)
    
    # # Regular case (no wrapping around)
    # return (starti <= endj) and (endi >= startj)

# Delete all genes found in blocks from the genomes
def deleteGenes(genome1, genome2, blocks):
    for i in range(len(blocks)):
        for j in range(blocks[i][0],blocks[i][1]+1):
            genome1[j%len(genome1)] = -1
        for j in range(min(blocks[i][2],blocks[i][3]),max(blocks[i][2],blocks[i][3])+1):
            genome2[j%len(genome2)] = -1

# Convert block count to frequency count
def normalizeMatrix(blockCount):
    # Delete useless entries
    while (blockCount[-1] == 0):
        blockCount.pop()

    # Calculate the total amount of blocks
    totalBlocks = 0
    for i in range(len(blockCount)):
        totalBlocks += blockCount[i]
    #print(totalBlocks)

    if totalBlocks > 0:
        for i in range(len(blockCount)):
            blockCount[i] /= totalBlocks

# Compare 2 distributions
def compareDistributions(distributionReal, distributionSimulation):
    # Sort the real distribution into bins as equally as possible (max 10 bins)
    bins = [0]
    bestBins = bins
    binMaxLength = [0] # The maximum block length (minus 1 because arrays start at 0) found in the bin i
    bestBinMaxLength = binMaxLength
    # Try 10 bins, then 9 bins, then 8... until 2 bins and find the number that gives the best ratio, give some leeway to bigger bin counts
    for binCount in range(10, 1, -1):
        # Reset the current array of bins
        bins = [0] * binCount
        binMaxLength = [-1] * binCount
        # Handle bins
        for i in range(0, binCount):
            for j in range(max(binMaxLength)+1, len(distributionReal)):
                bins[i] += distributionReal[j]
                binMaxLength[i] = j
                if(bins[i] >= 1/binCount):
                    break
            
        if (min(bestBins) == 0 or (min(bins) > 0 and max(bins)/min(bins) < (max(bestBins)/min(bestBins)) - 0.01)):
            bestBins = bins
            bestBinMaxLength = binMaxLength
    
    # Remove bins of capacity 0
    for i in range(len(bestBins)-1, 0, -1):
        if (bestBins[i] == 0):
            bestBins.pop(i)
            bestBinMaxLength.pop(i)

    # Sort the simulated distribution into bins following the ranges we got from the real data
    binsSimulation = [0] * len(bestBins)
    for j in range(0, bestBinMaxLength[0]+1):
        binsSimulation[0] += distributionSimulation[j]
    for i in range(1, len(bestBins)):
        for j in range(bestBinMaxLength[i-1]+1, bestBinMaxLength[i]+1):
            binsSimulation[i] += distributionSimulation[j]
        # Add all remaining blocks to the last bin
        if (i == len(bestBins)-1):
            for j in range(bestBinMaxLength[i]+1, len(distributionSimulation)):
                binsSimulation[i] += distributionSimulation[j]

    # print("Synteny block lengths bin distribution:")
    # print("1-" + str(bestBinMaxLength[0]+1), end= " ")
    # for i in range(1, len(bestBins)):
    #     print(str(bestBinMaxLength[i-1]+2) + "-" + str(bestBinMaxLength[i]+1), end= " ")
    # print("Real data bins:")
    # print(bestBins)
    # print("Simulation bins:")
    # print(binsSimulation)

    return wasserstein_distance(bestBins, binsSimulation) 

# Run the simulation many times and take the average
def runSimulation():
    GENOME_SIZE = int(jumpgenomesize.get())
    TIMES = int(times.get())
    totalBlockCount = [0] * GENOME_SIZE
    for i in range(TIMES):
        # Initiate the first genome
        genome1 = [0] * GENOME_SIZE
        for j in range(GENOME_SIZE):
            genome1[j] = j
        
        # Get the second genome using the Poisson jump model
        genome2 = jump(genome1, float(edgelength.get()))

        # Compare the 2 genomes
        blockCount = findSyntenyJump(genome1,genome2)
        for j in range(GENOME_SIZE):
            totalBlockCount[j] += blockCount[j]

    display_results(totalBlockCount)

# Compare 2 manually input genomes
def compareGenomes():
    genome1 = list(map(int, str(gen1.get()).split()))
    genome2 = list(map(int, str(gen2.get()).split()))

    # Clean singletons
    genomeCompare1 = [0] * (len(genome1) + len(genome2))
    genomeCompare2 = [0] * (len(genome1) + len(genome2))
    for i in range(max(len(genome1), len(genome2))):
        if (i < len(genome1)):
            genomeCompare1[genome1[i]-1] = 1
        if (i < len(genome2)):
            genomeCompare2[genome2[i]-1] = 1

    for i in range(len(genomeCompare1)):
        # Gene exists in genome1 but not in genome2
        if (genomeCompare1[i] > genomeCompare2[i]):
            for j in range(len(genome1)-1, -1, -1):
                if (genome1[j] == i+1):
                    genome1.pop(j)
        # Gene exists in genome2 but not in genome1
        elif (genomeCompare1[i] < genomeCompare2[i]):
            for j in range(len(genome2)-1, -1, -1):
                if (genome2[j] == i+1):
                    genome2.pop(j)

    print("Comparing genomes:\n" + str(genome1) + "\n" + str(genome2))
    blocks = findSyntenyReal(genome1,genome2)
    print(blocks)

# Extract genome data
def extractGeneFromRD(line):
    splitLine = line.split()
    genomeId= splitLine[3]
    genomeCog = int(splitLine[-1][4:])
    partitionId= splitLine[4] # NC_ number
    return genomeId, genomeCog, partitionId
 
# Read the file and extract genome data
def runFromFile(filepath):
    if not filepath:
        print("No file selected.")
        return
    
    genomes = []
    print(filepath)
    global FILE_NAME
    # FILE_NAME = filepath.split("/")[:8]
    FILE_NAME = os.path.splitext(os.path.basename(filepath))[0]
    
    with open(filepath, "r") as file:
        lines = file.readlines()
        genomeNames = []
        partitions = [] # format = (name, length)
        curr_id = ""
        curr_part = ""
        for line in lines:
            if line[0] == "#":
                continue
            new_id, new_num, new_part = extractGeneFromRD(line)
            if new_id != curr_id or curr_part != new_part:
                genomes.append([])
                genomeNames.append(new_id)
                partitions.append(new_part)
                curr_id = new_id
                curr_part = new_part
            genomes[-1].append(new_num)

        # Print the matrix
        # for genome in genomes:
        #     print(genome)
        return genomes, genomeNames, partitions

# Generate graph for the ATGC from read file
def runATGC(genomes, genomeNames, partitions):
    # genomes, genomeNames, partitions = runFromFile()
    if genomes is None:
        print("Error: Empty matrix")
        return
    
    # Make output folder
    os.makedirs("output", exist_ok=True)

    # Find max length
    genomeSize = 0
    for i in range(len(genomes)):
        if genomeSize < len(genomes[i]):
            genomeSize = len(genomes[i])

    # Compare every pair of genomes and write the results in a file
    with open(os.path.join("output", FILE_NAME + ".tab"), 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')

        # Group all partitions of each genome into 1 group
        genomeGroups = {}
        for idx, name in enumerate(genomeNames):
            if name not in genomeGroups:
                genomeGroups[name] = []
            genomeGroups[name].append(idx)
        
        uniqueNames = list(genomeGroups.keys())

        pairs_to_compare = []
        for i in range(len(uniqueNames)):
            for j in range(i+1, len(uniqueNames)):
                # Get all partitions for both genomes
                group1 = genomeGroups[uniqueNames[i]]
                group2 = genomeGroups[uniqueNames[j]]
                
                # Generate all combinations between partitions
                for pair in product(group1, group2):
                    pairs_to_compare.append(pair)

        # Remnant before we did grouping
        # for i in range(0, len(genomes)):
        #     for j in range(i+1, len(genomes)):
        #         if(genomeNames[i] != genomeNames[j]):
        #             pairs_to_compare.append((i, j))
        
        count = 1
        for i,j in pairs_to_compare:
            # detailed printing to test partition sort
            # print("Comparing ", genomeNames[i], "partition", partitions[i], "with ", genomeNames[j], "partition", partitions[j], "of length", len(genomes[i]) * len(genomes[j]))
            if(count%10 ==0):
                print(count,"out of",len(pairs_to_compare))
            count += 1
            # Compare genome i with genome j
            # print("\t", j-i, "out of", len(genomes)-i-1)
            blocks = findSyntenyReal(genomes[i].copy(),genomes[j].copy())

            # Write to tab file
            writer.writerow([genomeNames[i], genomeNames[j], partitions[i], partitions[j], len(genomes[i]), len(genomes[j])])
            writer.writerows(blocks)
            writer.writerow([])
    
    print("finished")

# Browse file
def browse_file():
    # Ask if the user wants to select a single file or an entire folder
    dialog = tk.Toplevel()
    dialog.title("Select File or Folder")
    dialog.geometry("400x150")
    
    def select_file():
        path = filedialog.askopenfilename(filetypes=(("Pseudo-terminal utilities", "*.pty"), ("All files", "*.*")))
        if path:
            file_entry.delete(0, tk.END)
            file_entry.insert(0, path)
        dialog.destroy()
    
    def select_folder():
        path = filedialog.askdirectory()
        if path:
            file_entry.delete(0, tk.END)
            file_entry.insert(0, path)
        dialog.destroy()
    
    ttk.Label(dialog, text="Choose selection type:").pack(pady=10)
    ttk.Button(dialog, text="📄 Select File", command=select_file).pack(fill='x', padx=50, pady=5)
    ttk.Button(dialog, text="📁 Select Folder", command=select_folder).pack(fill='x', padx=50, pady=5)

# Run the function on a single file or every file in a folder
def runFile():
    path = file_entry.get()
    if not path:
        print("No file or folder selected.")
        return
    
    # If a file is selected, process it
    if os.path.isfile(path):
        print(f"Processing: {os.path.basename(path)}")
        genomes, genomeNames, partitions = runFromFile(path)
        runATGC(genomes, genomeNames, partitions)
    
    # If a folder is selected, process all files in it
    elif os.path.isdir(path):
        for filename in sorted(os.listdir(path)):
            filepath = os.path.join(path, filename)
            if os.path.isfile(filepath):
                try:
                    print(f"Processing: {filename}")
                    genomes, genomeNames, partitions = runFromFile(filepath)
                    runATGC(genomes, genomeNames, partitions)
                except Exception as e:
                    print(f"Error processing {filename}: {str(e)}")
    
    else:
        print("Invalid path")


# GUI setup
root = tk.Tk()
root.title("Genome Synteny Simulation")
root.geometry("900x600")
frame_left = tk.Frame(root, padx=10, pady=10)
frame_left.pack(side=tk.LEFT, fill=tk.Y)
title_label = tk.Label(frame_left, text="Genome Synteny Simulator", font=("Arial", 16))
title_label.pack(pady=10)

# Edge length
ttk.Label(frame_left, text="Edge Length").pack()
edgelength = ttk.Entry(frame_left)
edgelength.insert(0, "0.2")  # Default value
edgelength.pack(pady=(0,30))

# Read file
ttk.Label(frame_left, text="ATGC Data File").pack()
file_entry = ttk.Entry(frame_left)
file_entry.pack()
browse_button = ttk.Button(frame_left, text="Browse", command=browse_file)
browse_button.pack()
file_sim_button = ttk.Button(frame_left, text="Find synteny blocks in file", command=runFile)
file_sim_button.pack(pady=(10,30))

# Jump simulation
ttk.Label(frame_left, text="Number of simulations").pack()
times = ttk.Entry(frame_left)
times.insert(0, "1")  # Default value
times.pack()
ttk.Label(frame_left, text="Genome size (jump sim)").pack()
jumpgenomesize = ttk.Entry(frame_left)
jumpgenomesize.insert(0, "5000")  # Default value
jumpgenomesize.pack()
run_button = ttk.Button(frame_left, text="Run Jump Simulation", command=runSimulation)
run_button.pack(pady=(10,30))

# Direct comparison of 2 genomes
ttk.Label(frame_left, text="Genomes").pack()
gen1 = ttk.Entry(frame_left)
gen1.insert(0, "")  # Default value
gen1.pack()
gen2 = ttk.Entry(frame_left)
gen2.insert(0, "")  # Default value
gen2.pack()
run_button = ttk.Button(frame_left, text="Direct comparison", command=compareGenomes)
run_button.pack(pady=(10,0))

# Note
color_guide = tk.Label(frame_left, text="Red = formula\nBlue = induced")
color_guide.pack(pady=10)

frame_right = tk.Frame(root, padx=10, pady=10)
frame_right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

root.mainloop()