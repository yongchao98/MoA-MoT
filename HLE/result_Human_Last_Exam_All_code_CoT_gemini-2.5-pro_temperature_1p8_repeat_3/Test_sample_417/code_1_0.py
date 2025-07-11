# This script calculates the expected number of signals in the 1H NMR spectrum
# of 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene
# based on its molecular symmetry and structure.

# The molecule possesses C3 rotational symmetry. This simplifies the problem by making
# the three large substituent "arms" and the three central methyl groups equivalent.
# We will count the unique proton environments on the central core and on a single arm.

# 1. Protons on the central 2,4,6-trimethylbenzene core.
# The three methyl groups are equivalent due to C3 symmetry.
signals_core_methyls = 1

# 2. Protons on one substituent "arm". The arm consists of a -CH2- linker
# and a chiral camphor-pyrazole group.

# 2a. Protons on the -CH2- linker.
# The CH2 group is adjacent to a chiral center, making its two protons diastereotopic
# and chemically non-equivalent.
signals_linker_ch2 = 2

# 2b. Protons on the rigid, chiral camphor-pyrazole moiety.
# We count each unique proton environment within this group.
signals_pyrazole_h = 1           # One proton on the pyrazole ring itself.
signals_bridgehead_h = 1         # One proton at the C4 bridgehead.
signals_c5_ch2 = 2               # Two non-equivalent protons (endo/exo) on the C5 methylene.
signals_c6_ch2 = 2               # Two non-equivalent protons (endo/exo) on the C6 methylene.
signals_c1_methyl = 1            # One methyl group at the C1 position.
signals_c7_gem_dimethyl = 2      # Two non-equivalent methyl groups at the C7 position.

# Summing up all the components to find the total number of unique signals.
all_components = [
    signals_core_methyls,
    signals_linker_ch2,
    signals_pyrazole_h,
    signals_bridgehead_h,
    signals_c5_ch2,
    signals_c6_ch2,
    signals_c1_methyl,
    signals_c7_gem_dimethyl
]

total_signals = sum(all_components)

# Printing the final breakdown and result.
print("Calculation of total expected 1H NMR peaks:")
print(f"Core Methyls: {signals_core_methyls}")
print(f"Linker -CH2-: {signals_linker_ch2}")
print(f"Pyrazole H: {signals_pyrazole_h}")
print(f"Bridgehead H: {signals_bridgehead_h}")
print(f"C5 -CH2- protons: {signals_c5_ch2}")
print(f"C6 -CH2- protons: {signals_c6_ch2}")
print(f"C1 -CH3 group: {signals_c1_methyl}")
print(f"C7 gem-dimethyl groups: {signals_c7_gem_dimethyl}")
print("-" * 30)

# Displaying the final equation with each individual count as requested.
equation_str = " + ".join(map(str, all_components))
print(f"Total Peaks = {equation_str}")
print(f"Total number of expected peaks = {total_signals}")
