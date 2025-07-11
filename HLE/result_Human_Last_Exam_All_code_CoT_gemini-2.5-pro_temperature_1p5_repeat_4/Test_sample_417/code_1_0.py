def count_nmr_peaks():
    """
    Calculates and explains the number of expected 1H NMR peaks for the given molecule.
    """
    print("Step-by-step analysis for determining the number of 1H NMR peaks:")

    # 1. Central Core Protons
    # The three methyl groups on the 2,4,6-trimethylbenzene core are equivalent due to C3 symmetry.
    core_peaks = 1
    print(f"\nNumber of peaks from the central core's three equivalent methyl groups: {core_peaks}")

    # 2. Methylene Linker Protons
    # The three -CH2- linkers are equivalent. Within each linker, the two protons
    # are diastereotopic because they are adjacent to a chiral center.
    linker_peaks = 2
    print(f"Number of peaks from the diastereotopic protons of the -CH2- linker: {linker_peaks}")

    # 3. Chiral Substituent Protons
    # The three chiral substituents are equivalent. We analyze one unit. It is asymmetric.
    # - 1 CH proton on the pyrazole ring
    # - 3 distinct methyl groups
    # - 1 CH proton at a bridgehead
    # - 2 CH2 groups, resulting in 4 non-equivalent protons
    substituent_ch_pyrazole = 1
    substituent_methyls = 3
    substituent_ch_bridgehead = 1
    substituent_ch2_protons = 4
    substituent_total_peaks = substituent_ch_pyrazole + substituent_methyls + substituent_ch_bridgehead + substituent_ch2_protons
    print(f"Number of peaks from one chiral substituent: {substituent_total_peaks}")
    print(f"  - Pyrazole CH: {substituent_ch_pyrazole}")
    print(f"  - Methyl groups: {substituent_methyls}")
    print(f"  - Bridgehead CH: {substituent_ch_bridgehead}")
    print(f"  - Bicyclic CH2 protons: {substituent_ch2_protons}")

    # 4. Total Calculation
    total_peaks = core_peaks + linker_peaks + substituent_total_peaks
    print("\nTotal number of peaks is the sum of all non-equivalent protons.")
    print("Final Equation:")
    print(f"{core_peaks} (core) + {linker_peaks} (linker) + {substituent_total_peaks} (substituent) = {total_peaks}")

    print(f"\nThus, there are {total_peaks} expected peaks in the 1H NMR spectrum.")

if __name__ == "__main__":
    count_nmr_peaks()
<<<H>>>