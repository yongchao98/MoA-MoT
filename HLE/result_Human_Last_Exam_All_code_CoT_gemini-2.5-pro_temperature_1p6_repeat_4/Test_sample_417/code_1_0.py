def count_nmr_peaks():
    """
    Calculates the number of expected 1H NMR peaks for the given molecule based on its structure and symmetry.
    """

    # Part 1: Protons on the central core and linkers
    # The three methyl groups at positions 2,4,6 are equivalent due to C3 symmetry.
    core_methyl_peaks = 1
    
    # The three methylene linkers at 1,3,5 are equivalent by C3 symmetry.
    # However, the two protons on each CH2 group are diastereotopic because they are attached to a chiral group.
    linker_methylene_peaks = 2

    # Part 2: Protons on one substituent arm
    # The arm is chiral and rigid, so we count non-equivalent protons.
    # (4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl
    
    # C3-H on pyrazole ring
    pyrazole_ch_peaks = 1
    
    # C4-H at a bridgehead
    bridgehead_ch_peaks = 1
    
    # C5-H2 and C6-H2 are two diastereotopic CH2 groups. All 4 protons are unique.
    bicyclic_methylene_peaks = 4
    
    # Methyl group at C7
    c7_methyl_peaks = 1
    
    # Two diastereotopic methyl groups at C8 (gem-dimethyl)
    c8_gem_dimethyl_peaks = 2
    
    # Sum of peaks from one arm
    arm_peaks = pyrazole_ch_peaks + bridgehead_ch_peaks + bicyclic_methylene_peaks + c7_methyl_peaks + c8_gem_dimethyl_peaks

    # Part 3: Total number of peaks
    total_peaks = core_methyl_peaks + linker_methylene_peaks + arm_peaks
    
    print("Determining the total number of 1H NMR peaks by summing the signals from chemically distinct protons:")
    print("\n1. Central Core and Linkers:")
    print(f"   - Signals from three equivalent Ar-CH3 groups: {core_methyl_peaks}")
    print(f"   - Signals from three equivalent Ar-CH2-N linkers (diastereotopic protons): {linker_methylene_peaks}")
    
    print("\n2. One Substituent Arm:")
    print(f"   - Sum of signals from one arm: {pyrazole_ch_peaks} (C3-H) + {bridgehead_ch_peaks} (C4-H) + {bicyclic_methylene_peaks} (C5/6-H2) + {c7_methyl_peaks} (C7-Me) + {c8_gem_dimethyl_peaks} (C8-Me2) = {arm_peaks}")

    print("\n3. Final Calculation:")
    print(f"   Total Peaks = (Core Methyls) + (Linker Methylenes) + (Arm Signals)")
    print(f"   Final Equation: {core_methyl_peaks} + {linker_methylene_peaks} + {arm_peaks} = {total_peaks}")

count_nmr_peaks()