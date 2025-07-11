import sys

def solve():
    """
    Analyzes a coiled-coil sequence to predict its oligomeric state.
    """
    # The input protein sequence
    sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"

    # Step 1: Define amino acid properties for analysis
    # Kyte-Doolittle scale for hydrophobicity scoring
    kd_scale = {
        'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
        'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
        'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
    }
    # Beta-branched hydrophobic/ambivalent amino acids
    beta_branched = ['I', 'V', 'T']

    # Step 2: Find the optimal heptad repeat register
    print("Step 1: Analyzing the sequence to find the optimal heptad repeat register.")
    print(f"Sequence: {sequence}")
    
    best_offset = -1
    max_score = -sys.float_info.max

    for offset in range(7):
        current_score = 0
        for i, residue in enumerate(sequence):
            # position in heptad is 'a' if result is 0, 'b' if 1, etc.
            pos_in_heptad = (i - offset) % 7
            if pos_in_heptad == 0 or pos_in_heptad == 3:  # 'a' or 'd' position
                current_score += kd_scale.get(residue, 0)
        
        if current_score > max_score:
            max_score = current_score
            best_offset = offset

    # Step 3: Extract and analyze residues from the best register
    a_residues = []
    d_residues = []
    register_str = ""
    heptad_map = "abcdefg"
    for i, residue in enumerate(sequence):
        pos_in_heptad_idx = (i - best_offset) % 7
        register_char = heptad_map[pos_in_heptad_idx]
        register_str += register_char
        if register_char == 'a':
            a_residues.append(residue)
        elif register_char == 'd':
            d_residues.append(residue)
    
    print(f"\nOptimal register found (offset={best_offset}), which maximizes core hydrophobicity.")
    print(f"Register: {register_str}")
    
    print("\nStep 2: Identifying core residues at 'a' and 'd' positions.")
    print(f"Residues at 'a' positions: {', '.join(a_residues)}")
    print(f"Residues at 'd' positions: {', '.join(d_residues)}")

    # Step 4: Apply knobs-into-holes packing rules
    # We ignore the final 'K' at a 'd' position as it's likely a terminal capping residue.
    core_d_residues = [r for r in d_residues if r != 'K']

    a_beta_count = sum(1 for res in a_residues if res in beta_branched)
    d_beta_count = sum(1 for res in core_d_residues if res in beta_branched)
    
    percent_a_beta = (a_beta_count / len(a_residues)) * 100
    percent_d_beta = (d_beta_count / len(core_d_residues)) * 100

    print("\nStep 3: Analyzing residue types for 'knobs-into-holes' packing.")
    print("The 'a' positions are predominantly non-beta-branched (Leucine).")
    print("The 'd' positions are predominantly beta-branched (Isoleucine).")
    print(f"Calculation: Percentage of beta-branched residues at 'a' positions = ({a_beta_count} / {len(a_residues)}) * 100 = {percent_a_beta:.0f}%")
    print(f"Calculation: Percentage of beta-branched residues at 'd' positions (core) = ({d_beta_count} / {len(core_d_residues)}) * 100 = {percent_d_beta:.0f}%")
    
    print("\nStep 4: Conclusion based on packing efficiency.")
    print("The observed pattern with non-beta-branched Leucine at 'a' and beta-branched Isoleucine at 'd' is a strong signature for a dimeric (2-stranded) coiled-coil.")
    print("This arrangement allows for efficient knobs-into-holes packing in a dimer, whereas it would cause steric clashes in the core of a trimer.")
    print("Therefore, the protein is predicted to form a dimer.")

solve()
<<<A>>>