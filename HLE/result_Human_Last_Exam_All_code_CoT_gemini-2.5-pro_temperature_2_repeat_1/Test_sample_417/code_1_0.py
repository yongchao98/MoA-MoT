def calculate_nmr_peaks():
    """
    This function analyzes the structure of 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-
    tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene and
    calculates the expected number of signals in its 1H NMR spectrum.
    """
    print("Analyzing the molecule to determine the number of 1H NMR peaks...")
    print("="*60)

    # 1. Central Core: 2,4,6-trimethylbenzene
    # The three methyl groups are equivalent due to C3 symmetry.
    central_core_signals = 1
    print(f"Contribution from central core's three equivalent methyl groups: {central_core_signals} peak")

    # 2. Linker: -CH2- groups
    # The two protons on the CH2 linker are diastereotopic due to the chiral substituent.
    linker_signals = 2
    print(f"Contribution from the diastereotopic protons of the -CH2- linker: {linker_signals} peaks")

    # 3. Indazole substituent
    # We count unique protons in one of the three equivalent substituent arms.
    # It's a rigid, asymmetric chiral system.
    
    # Proton on the pyrazole ring
    indazole_pyrazole_h = 1
    
    # Protons on the bicyclic (camphor-like) skeleton
    # Bridgehead methine proton (C4-H)
    camphor_methine_h = 1
    
    # Three methyl groups (C7-Me, and two at C8) are all non-equivalent
    camphor_methyls = 3
    
    # Four methylene protons (H-5exo, H-5endo, H-6exo, H-6endo) are non-equivalent
    camphor_methylenes = 4
    
    substituent_signals = indazole_pyrazole_h + camphor_methine_h + camphor_methyls + camphor_methylenes
    print("Contribution from one of the three equivalent indazole substituents:")
    print(f"  - Pyrazole ring H: {indazole_pyrazole_h} peak")
    print(f"  - Bridgehead methine H: {camphor_methine_h} peak")
    print(f"  - Three non-equivalent methyl groups: {camphor_methyls} peaks")
    print(f"  - Four non-equivalent methylene protons: {camphor_methylenes} peaks")
    print(f"  Total from substituent: {substituent_signals} peaks")

    print("="*60)
    
    # Total Calculation
    total_signals = central_core_signals + linker_signals + substituent_signals
    
    print("Final Calculation:")
    print(f"{central_core_signals} (from core methyls) + {linker_signals} (from CH2 linkers) + {substituent_signals} (from indazole unit) = {total_signals} peaks")
    print(f"\nThe total expected number of peaks in the 1H NMR spectrum is {total_signals}.")

calculate_nmr_peaks()