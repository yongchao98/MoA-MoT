def calculate_nmr_peaks():
    """
    Calculates the number of expected 1H NMR peaks for the given molecule based on its symmetry.
    """
    # 1. Central 2,4,6-trimethylbenzene ring
    # The three methyl groups are equivalent due to C3 symmetry.
    central_ring_methyl_signals = 1
    print(f"Number of signals from central ring methyls: {central_ring_methyl_signals}")

    # 2. Methylene (-CH2-) linker groups
    # The three -CH2- groups are equivalent due to C3 symmetry.
    # The two protons on each -CH2- are diastereotopic because they are attached to a chiral center.
    linker_ch2_signals = 2
    print(f"Number of signals from linker -CH2- groups: {linker_ch2_signals}")

    # 3. Chiral substituent arm
    # We only need to count signals for one arm, as the three arms are equivalent.
    # The arm is rigid and asymmetric, so all its protons are unique.
    h3_proton = 1
    h4_bridgehead_proton = 1
    c5_ch2_protons = 2  # Diastereotopic
    c6_ch2_protons = 2  # Diastereotopic
    non_equivalent_methyls = 3
    
    arm_signals = h3_proton + h4_bridgehead_proton + c5_ch2_protons + c6_ch2_protons + non_equivalent_methyls
    print(f"Number of signals from one chiral arm: {h3_proton} + {h4_bridgehead_proton} + {c5_ch2_protons} + {c6_ch2_protons} + {non_equivalent_methyls} = {arm_signals}")

    # 4. Total number of signals
    total_signals = central_ring_methyl_signals + linker_ch2_signals + arm_signals
    
    print("\n--- Total Calculation ---")
    print(f"Total number of peaks = {central_ring_methyl_signals} (ring methyls) + {linker_ch2_signals} (linker protons) + {arm_signals} (chiral arm)")
    print(f"Total number of peaks = {total_signals}")

calculate_nmr_peaks()