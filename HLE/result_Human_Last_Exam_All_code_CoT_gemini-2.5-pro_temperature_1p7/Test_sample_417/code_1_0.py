import collections

def calculate_nmr_peaks():
    """
    This script calculates the expected number of signals in the 1H NMR spectrum of
    1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene
    by analyzing its structure for chemically non-equivalent protons.
    """
    print("Determining the number of 1H NMR signals based on molecular structure:")

    # Using an ordered dictionary to maintain the order of analysis
    signal_sources = collections.OrderedDict()

    # 1. Analyze the central 2,4,6-trimethylbenzene core
    # The C3 symmetry makes the three methyl groups chemically equivalent.
    signal_sources["Central Ring -CH3 Groups"] = 1
    print("- The 3 methyl groups on the central benzene ring are equivalent due to C3 symmetry -> 1 signal.")

    # 2. Analyze one of the three identical substituent arms.
    # The arm is chiral, making diastereotopic protons/groups non-equivalent.
    print("- The 3 large substituents are identical. We analyze one arm:")

    # Benzylic -CH2- linker group
    signal_sources["Benzylic -CH2- Protons"] = 2
    print("  - The 2 protons of the benzylic -CH2- linker are diastereotopic -> 2 signals.")

    # Protons on the indazole ring system
    signal_sources["Pyrazole C3-H Proton"] = 1
    print("  - The proton on the C3 position of the pyrazole ring -> 1 signal.")

    # Protons on the bicyclic framework
    signal_sources["Bridgehead C4-H Proton"] = 1
    print("  - The bridgehead proton at C4 -> 1 signal.")
    
    signal_sources["C5 Methylene Protons"] = 2
    print("  - The 2 protons at C5 are diastereotopic (endo/exo) -> 2 signals.")
    
    signal_sources["C6 Methylene Protons"] = 2
    print("  - The 2 protons at C6 are also diastereotopic (endo/exo) -> 2 signals.")

    # Methyl groups within the arm
    signal_sources["C7-Methyl Group"] = 1
    print("  - The methyl group at the C7 position -> 1 signal.")

    signal_sources["C8-gem-dimethyl Groups"] = 2
    print("  - The 2 methyl groups at C8 (gem-dimethyl) are diastereotopic -> 2 signals.")

    # 3. Sum up all the signals
    total_signals = sum(signal_sources.values())
    
    print("\nTotal number of signals is the sum of signals from all unique environments:")
    
    equation_parts = [str(v) for v in signal_sources.values()]
    equation_string = " + ".join(equation_parts)
    
    print(f"Final Equation: {equation_string} = {total_signals}")

calculate_nmr_peaks()