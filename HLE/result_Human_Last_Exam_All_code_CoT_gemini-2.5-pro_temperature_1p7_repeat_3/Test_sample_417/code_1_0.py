#
# Plan to calculate the number of 1H NMR signals for the given molecule.
# 1. Identify the overall molecular symmetry to find equivalent parts of the molecule.
# 2. Count the number of unique proton signals on the central benzene core.
# 3. Count the number of unique proton signals in one of the substituent arms, as they are all equivalent.
# 4. Sum the counts to find the total number of expected peaks.
#

# Step 1: Account for symmetry. The molecule has C3 symmetry.

# Step 2: Count signals from the central 2,4,6-trimethylbenzene core.
# The three methyl groups are equivalent due to C3 symmetry.
num_ar_me_signals = 1

# Step 3: Count signals from one substituent arm.
# The arm is a -CH2- linker attached to a chiral camphor-indazole ligand.

# The two protons of the -CH2- linker are diastereotopic due to the adjacent chiral ligand.
num_linker_ch2_signals = 2

# The camphor-indazole ligand is chiral and asymmetric. We count its unique protons.
# a) One proton on the pyrazole ring.
num_pyrazole_h_signals = 1

# b) Three methyl groups on the camphor framework become non-equivalent due to loss of symmetry.
num_camphor_me_signals = 3

# c) One proton at the bridgehead of the camphor framework.
num_camphor_bridgehead_ch_signals = 1

# d) Two CH2 groups on the camphor framework contain four non-equivalent protons in total.
num_camphor_methylene_ch_signals = 4

# Step 4: Sum all unique signals.
total_signals = (num_ar_me_signals +
                 num_linker_ch2_signals +
                 num_pyrazole_h_signals +
                 num_camphor_me_signals +
                 num_camphor_bridgehead_ch_signals +
                 num_camphor_methylene_ch_signals)

# Print the final equation with each component's contribution.
print("The total number of expected peaks is the sum of signals from each unique proton environment.")
print(f"Total peaks = (Ar-Me) + (Linker CH2) + (Pyrazole H) + (Camphor Me) + (Camphor CH) + (Camphor CH2s)")
print(f"Total peaks = {num_ar_me_signals} + {num_linker_ch2_signals} + {num_pyrazole_h_signals} + {num_camphor_me_signals} + {num_camphor_bridgehead_ch_signals} + {num_camphor_methylene_ch_signals}")
print(f"Total expected peaks = {total_signals}")