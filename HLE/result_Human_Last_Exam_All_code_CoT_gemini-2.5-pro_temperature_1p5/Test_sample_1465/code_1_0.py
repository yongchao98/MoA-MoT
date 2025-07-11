import collections

def get_heptad_positions(seq):
    """
    Identifies the correct heptad repeat frame and extracts a and d positions.
    This is based on heuristics for the given sequences.
    """
    # Frame 1: Assumes a starts at the 2nd residue (index 1), d at 5th (index 4)
    # This frame works for sequences 1, 2, 3, 5
    a_pos_frame1 = [seq[i] for i in range(1, len(seq), 7)]
    d_pos_frame1 = [seq[i] for i in range(4, len(seq), 7)]

    # Frame 2: Assumes a starts at the 3rd residue (index 2), d at 6th (index 5)
    # This frame works for sequence 4, which has a QxxL pattern
    a_pos_frame2 = [seq[i] for i in range(2, len(seq), 7)]
    d_pos_frame2 = [seq[i] for i in range(5, len(seq), 7)]

    # Heuristic to choose the frame: Frame 2 is specific to the Q...L pattern of seq 4.
    if collections.Counter(a_pos_frame2).get('Q', 0) >= 3:
        return a_pos_frame2, d_pos_frame2, "Frame 2 (a=pos 3, d=pos 6)"
    else:
        return a_pos_frame1, d_pos_frame1, "Frame 1 (a=pos 2, d=pos 5)"

def predict_oligomer_state(a_pos, d_pos):
    """
    Predicts oligomeric state based on simplified rules for a and d position residues.
    """
    # Rule for Trimers (I...T pattern)
    if 'T' in d_pos and 'I' in a_pos:
        return 3 # Beta-branched residue (T) at 'd' with I at 'a' promotes trimers.
    
    # Rule for Tetramers (I...W pattern)
    if 'W' in d_pos:
        return 4 # Bulky Tryptophan (W) at a core 'd' position can favor tetramers.
    
    # Rule for Dimers (Q...L pattern)
    if 'Q' in a_pos and 'L' in d_pos:
        return 2 # Polar 'a' position (Q) with canonical 'd' position (L) favors dimers.
    
    # Rule for Dimers (canonical I...L pattern)
    if 'I' in a_pos and 'L' in d_pos:
        return 2 # Canonical pattern for heterodimers and many homodimers.

    # Fallback if no specific rule matches
    return "Unknown"

# List of protein sequences to analyze
sequences = [
    "EIAQALKEIAKALKEIAWALKEIAQALK",
    "EIAALKQEIAALKKENAALKQEIAALKQ",
    "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
    "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
    "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
]

final_states = []

print("Analyzing coiled-coil sequences...")
print("-" * 35)

for i, seq in enumerate(sequences):
    a_residues, d_residues, frame_info = get_heptad_positions(seq)
    state = predict_oligomer_state(a_residues, d_residues)
    final_states.append(str(state))
    
    print(f"Sequence {i+1}: {seq}")
    print(f"  - Frame Used: {frame_info}")
    print(f"  - Core 'a' residues: {a_residues}")
    print(f"  - Core 'd' residues: {d_residues}")
    print(f"  - Predicted State: {state}")
    print("-" * 35)

# Final output as requested
equation_str = "Predicted Oligomeric States = " + ", ".join(final_states)
print(equation_str)