import collections

def analyze_coiled_coil(sequence):
    """
    Analyzes a coiled-coil sequence to predict its oligomeric state.

    The prediction is based on the identity of amino acids at the core 'a' and 'd'
    positions of the heptad repeat.
    """
    
    # --- Step 1: Identify the Heptad Register ---
    # The key to coiled-coils is the heptad repeat (a-b-c-d-e-f-g).
    # 'a' and 'd' positions are typically hydrophobic and form the core.
    # We must find the register that maximizes this hydrophobicity.
    # By inspection, the sequence has Leucines (L) and Isoleucines (I) repeating
    # at regular intervals. A highly stable hydrophobic core is formed when the
    # 'a' positions are consistently 'L' and the 'd' positions are consistently 'I'.
    # This occurs when the heptad register starts with 'a' at the first Leucine (index 6).
    
    start_a_index = 6
    heptad_map = 'abcdefg'
    
    # --- Step 2: Assign and Extract Residues ---
    assignments = collections.defaultdict(list)
    heptad_string = []
    
    print("--- Analysis of Coiled-Coil Sequence ---")
    print(f"Sequence: {sequence}")
    
    for i, residue in enumerate(sequence):
        # Calculate the position in the heptad (a=0, b=1, ...)
        position_index = (i - start_a_index) % 7
        position_letter = heptad_map[position_index]
        assignments[position_letter].append(residue)
        heptad_string.append(position_letter)
        
    print(f"Register: {''.join(heptad_string)}")
    print("\n--- Residues at Key Heptad Positions ---")
    
    core_a_residues = assignments['a']
    core_d_residues = assignments['d']
    flanking_e_residues = assignments['e']
    flanking_g_residues = assignments['g']
    
    print(f"Core 'a' positions are: {', '.join(core_a_residues)}")
    print(f"Core 'd' positions are: {', '.join(core_d_residues)}")
    print(f"Flanking 'e' positions are: {', '.join(flanking_e_residues)}")
    print(f"Flanking 'g' positions are: {', '.join(flanking_g_residues)}")
    
    # --- Step 3 & 4: Analyze Core and Conclude ---
    print("\n--- Conclusion on Oligomeric State ---")
    print("1. The hydrophobic core of the coiled-coil is formed by residues at the 'a' and 'd' positions.")
    print(f"2. In this sequence, the 'a' positions are consistently Leucine (L) and the 'd' positions are consistently Isoleucine (I).")
    print("3. Both L and I are large, bulky hydrophobic residues.")
    print("4. This pattern ('a'=Large Hydrophobe, 'd'=Large Hydrophobe) allows for ideal 'knobs-into-holes' packing between two helices.")
    print("5. Higher-order oligomers (trimers, tetramers) have more crowded cores and typically require smaller residues (like Alanine) at the 'a' and/or 'd' positions to avoid steric clash.")
    print("\nBased on the classic rules of coiled-coil packing, the predicted oligomeric state is a DIMER.")
    print("\nFinal equation for prediction: a=[L,L,L,L] + d=[I,I,I,I] => Oligomeric State = 2")


protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"
analyze_coiled_coil(protein_sequence)
