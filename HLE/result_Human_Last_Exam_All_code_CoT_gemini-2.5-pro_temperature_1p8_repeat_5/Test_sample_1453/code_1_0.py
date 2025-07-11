import collections

def analyze_coiled_coil(sequence):
    """
    Analyzes a protein sequence to determine its likely coiled-coil oligomeric state.
    """
    print(f"Analyzing sequence: {sequence}\n")

    # Step 1: Find the best heptad repeat register
    hydrophobic_residues = ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
    best_register = -1
    max_score = -1
    register_data = {}

    for start_index in range(7):
        a_pos_residues = []
        d_pos_residues = []
        current_score = 0
        
        # 'a' positions are at index i where (i - start_index) % 7 == 0
        # 'd' positions are at index i where (i - start_index) % 7 == 3
        
        for i, residue in enumerate(sequence):
            if (i - start_index) % 7 == 0:
                a_pos_residues.append(residue)
                if residue in hydrophobic_residues:
                    current_score += 1
            elif (i - start_index) % 7 == 3:
                d_pos_residues.append(residue)
                if residue in hydrophobic_residues:
                    current_score += 1

        register_data[start_index] = {
            'score': current_score,
            'a': a_pos_residues,
            'd': d_pos_residues
        }
        
        if current_score > max_score:
            max_score = current_score
            best_register = start_index

    print("--- Step 1 & 2: Identify Register and Extract Core Residues ---")
    print(f"The most likely heptad repeat starts at position {best_register + 1} of the sequence.")
    
    a_residues = register_data[best_register]['a']
    d_residues = register_data[best_register]['d']
    
    print(f"Residues at the 'a' positions: {a_residues}")
    print(f"Residues at the 'd' positions: {d_residues}\n")

    # Step 3 & 4: Apply Packing Principles and Conclude
    print("--- Step 3 & 4: Apply Packing Principles to Determine Oligomeric State ---")
    
    a_counts = collections.Counter(a_residues)
    d_counts = collections.Counter(d_residues)

    # Simplified rule-based deduction for common patterns
    conclusion = ""
    oligomeric_state = "Unknown"

    # Check for the pattern a=Ala, d=Leu
    is_a_Ala = all(res == 'A' for res in a_residues)
    is_d_Leu = all(res == 'L' for res in d_residues)

    if is_a_Ala and is_d_Leu:
        oligomeric_state = 2
        conclusion = (
            "The hydrophobic core is composed of Alanine (A) at the 'a' positions and Leucine (L) at the 'd' positions.\n"
            "This 'a=A, d=L' pattern is a classic hallmark of a dimeric coiled-coil.\n"
            "Reasoning (Knobs-into-Holes): The small side chain of Alanine at the 'a' position ('knob') packs efficiently into the 'hole' of a two-helix bundle (dimer).\n"
            "In a trimer or tetramer, the central core would be too large, and the small Alanine side chains would lead to poor, unstable packing."
        )
    else:
        conclusion = "The residue pattern does not match a simple, canonical coiled-coil rule. Further analysis would be required."

    print(conclusion)
    
    if oligomeric_state != "Unknown":
        print("\nFinal Determination:")
        print(f"Oligomeric State = {oligomeric_state}")

# The sequence in question
protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"
analyze_coiled_coil(protein_sequence)