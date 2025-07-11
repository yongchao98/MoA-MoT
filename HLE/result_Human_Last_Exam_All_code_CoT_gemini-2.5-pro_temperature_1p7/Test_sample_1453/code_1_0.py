import collections

def analyze_coiled_coil_oligomerization():
    """
    Analyzes a protein sequence to determine its likely coiled-coil oligomeric state
    by identifying the heptad repeat and analyzing the core 'a' and 'd' positions.
    """
    sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"

    # A simple set of common hydrophobic core residues.
    hydrophobic_residues = {'L', 'I', 'V', 'A', 'M', 'F', 'W'}

    best_frame_start_index = -1
    max_hydrophobic_score = -1
    best_a_residues = []
    best_d_residues = []

    # 1. Test all 7 possible frames to find the optimal heptad repeat.
    for frame_start in range(7):
        a_positions = []
        d_positions = []
        hydrophobic_score = 0
        
        for i in range(len(sequence)):
            # Determine the position (0-6) in the heptad repeat (a-g)
            heptad_pos = (i - frame_start + 7) % 7
            
            # Position 'a' corresponds to index 0, 'd' to index 3
            if heptad_pos == 0:
                a_positions.append(sequence[i])
                if sequence[i] in hydrophobic_residues:
                    hydrophobic_score += 1
            elif heptad_pos == 3:
                d_positions.append(sequence[i])
                if sequence[i] in hydrophobic_residues:
                    hydrophobic_score += 1

        # The best frame is the one with the most hydrophobic residues in the core
        if hydrophobic_score > max_hydrophobic_score:
            max_hydrophobic_score = hydrophobic_score
            best_frame_start_index = frame_start
            best_a_residues = a_positions
            best_d_residues = d_positions

    print("--- Coiled-Coil Analysis ---")
    print(f"Sequence: {sequence}")
    print("\nStep 1: The optimal heptad repeat was identified.")
    
    # 2. Display the residues at the identified core positions.
    print("\nStep 2: The core residues for the optimal repeat are:")
    # We must output each residue in the 'equation' as requested
    print("Residues at 'a' positions: " + ", ".join(best_a_residues))
    print("Residues at 'd' positions: " + ", ".join(best_d_residues))

    # 3. Apply biochemical rules to predict the oligomeric state.
    print("\nStep 3: Predicting the oligomeric state based on core composition.")
    
    # Use Counter to find the most common residue at each position
    dominant_a = collections.Counter(best_a_residues).most_common(1)[0][0]
    dominant_d = collections.Counter(best_d_residues).most_common(1)[0][0]

    final_state = 0
    explanation = ""

    if dominant_a == 'A' and dominant_d == 'L':
        final_state = 2
        explanation = "This composition (Alanine at 'a', Leucine at 'd') is a canonical motif for a Dimer (2-stranded coil).\nThe small side chain of Alanine packs efficiently against the larger Leucine side chain, creating a stable hydrophobic core for a dimer."
    elif dominant_a in ('I', 'L') and dominant_d in ('I', 'L'):
        final_state = 3
        explanation = "A core composed primarily of bulky Isoleucine or Leucine residues often favors a Trimer (3-stranded coil) for optimal packing."
    else:
        final_state = "Undetermined by simple rules"
        explanation = "The residue pattern does not fit a common oligomerization motif."

    # Final "equation" showing the logic and result
    print(f"\nFinal Equation: a=[{', '.join(best_a_residues)}] + d=[{', '.join(best_d_residues)}] => Oligomeric State = {final_state}")
    print("\nExplanation: " + explanation)

# Execute the analysis
analyze_coiled_coil_oligomerization()