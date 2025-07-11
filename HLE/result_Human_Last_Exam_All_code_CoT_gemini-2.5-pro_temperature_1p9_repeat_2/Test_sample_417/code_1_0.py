def count_nmr_peaks():
    """
    Calculates the number of expected 1H NMR peaks for the given molecule based on its structure and symmetry.
    """
    # 1. Signals from the central mesitylene core.
    # The three methyl groups are equivalent due to C3 symmetry.
    num_signals_core = 1
    
    # 2. Signals from the -CH2- linking group in one arm.
    # The two protons are diastereotopic because the attached ligand is chiral.
    num_signals_ch2_bridge = 2
    
    # 3. Signals from one of the chiral ligands.
    # The ligand is a complex, rigid bicyclic system.
    
    # Proton on the pyrazole ring (C3-H)
    signals_pyrazole_h = 1
    
    # Three non-equivalent methyl groups (one on C7, two diastereotopic on C8)
    signals_ligand_methyls = 3
    
    # Protons on the saturated bicyclic frame:
    # Bridgehead proton at C4
    signals_bridgehead_h = 1
    # Two diastereotopic protons on C5
    signals_c5_h = 2
    # Two diastereotopic protons on C6
    signals_c6_h = 2
    
    num_signals_ligand = signals_pyrazole_h + signals_ligand_methyls + signals_bridgehead_h + signals_c5_h + signals_c6_h
    
    # 4. Total number of signals is the sum from all unique parts.
    total_signals = num_signals_core + num_signals_ch2_bridge + num_signals_ligand

    print("Analysis of 1H NMR signals:")
    print("-" * 30)
    print(f"Signals from central core's equivalent methyl groups: {num_signals_core}")
    print(f"Signals from one arm's diastereotopic -CH2- bridge protons: {num_signals_ch2_bridge}")
    print(f"Signals from one arm's chiral ligand protons: {num_signals_ligand}")
    print("-" * 30)
    print("The final equation is the sum of signals from each unique part:")
    print(f"Total Signals = {num_signals_core} (core) + {num_signals_ch2_bridge} (bridge) + {num_signals_ligand} (ligand) = {total_signals}")

count_nmr_peaks()
