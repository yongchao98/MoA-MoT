def count_nmr_signals():
    """
    Calculates the number of expected 1H NMR peaks for the given molecule
    by analyzing its structure and symmetry.
    """
    
    # Analysis based on molecular structure and symmetry (C3)
    
    # 1. Signals from the central 2,4,6-trimethylbenzene core.
    # The three methyl groups are equivalent due to C3 symmetry.
    signals_core_methyls = 1
    
    # 2. Signals from the three equivalent methylene (-CH2-) bridges.
    # The two protons on each bridge are diastereotopic due to the adjacent
    # chiral group, making them non-equivalent.
    signals_bridge_ch2 = 2
    
    # 3. Signals from one of the three equivalent chiral camphopyrazole arms.
    # We only need to count the unique protons within one arm.
    
    # Pyrazole ring proton
    signals_pyrazole_h = 1
    
    # Protons on the rigid camphor-like skeleton
    # Tertiary bridgehead proton (-CH-)
    signals_camphor_ch = 1
    
    # Two different methylene groups (-CH2-). Each has 2 non-equivalent protons.
    signals_camphor_ch2_a = 2
    signals_camphor_ch2_b = 2
    
    # Three methyl groups (-CH3). All are chemically distinct due to the
    # asymmetric environment.
    signals_camphor_methyls = 3
    
    # Calculate the total number of signals by summing the components.
    total_signals = (
        signals_core_methyls +
        signals_bridge_ch2 +
        signals_pyrazole_h +
        signals_camphor_ch +
        signals_camphor_ch2_a +
        signals_camphor_ch2_b +
        signals_camphor_methyls
    )

    print("Step-by-step count of 1H NMR signals:")
    print(f"1. Central core methyls (equivalent by symmetry): {signals_core_methyls}")
    print(f"2. Methylene bridge protons (diastereotopic): {signals_bridge_ch2}")
    print(f"3. Pyrazole ring proton: {signals_pyrazole_h}")
    print(f"4. Camphor skeleton -CH- proton: {signals_camphor_ch}")
    print(f"5. Camphor skeleton -CH2- group 'A' (diastereotopic): {signals_camphor_ch2_a}")
    print(f"6. Camphor skeleton -CH2- group 'B' (diastereotopic): {signals_camphor_ch2_b}")
    print(f"7. Camphor skeleton methyl groups (all non-equivalent): {signals_camphor_methyls}")
    print("\nFinal calculation:")
    print(f"Total signals = {signals_core_methyls} + {signals_bridge_ch2} + {signals_pyrazole_h} + {signals_camphor_ch} + {signals_camphor_ch2_a} + {signals_camphor_ch2_b} + {signals_camphor_methyls} = {total_signals}")

count_nmr_signals()