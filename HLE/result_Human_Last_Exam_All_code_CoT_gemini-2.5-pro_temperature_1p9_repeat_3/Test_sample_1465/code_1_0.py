import collections

def solve():
    """
    Predicts the oligomeric state of coiled-coil sequences based on their
    hydrophobic core composition.
    """
    sequences = [
        "EIAQALKEIAKALKEIAWALKEIAQALK",
        "EIAALKQEIAALKKENAALKQEIAALKQ",
        "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
        "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
        "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
    ]

    # Simplified hydrophobicity scale for phasing
    # (higher value = more hydrophobic)
    hydrophobicity = collections.defaultdict(int, {
        'I': 3, 'L': 3, 'V': 3, 'W': 3, 'F': 3,
        'A': 2, 'M': 2, 'T': 1
    })

    def find_best_frame(sequence):
        """
        Finds the best heptad repeat frame by maximizing the hydrophobicity
        of the core 'a' and 'd' positions.
        """
        best_frame = None
        max_score = -1

        for frame_start in range(7):
            a_pos_residues = []
            d_pos_residues = []
            current_score = 0

            for i, residue in enumerate(sequence):
                pos = (i - frame_start) % 7
                if pos == 0:  # 'a' position
                    a_pos_residues.append(residue)
                    current_score += hydrophobicity[residue]
                elif pos == 3:  # 'd' position
                    d_pos_residues.append(residue)
                    current_score += hydrophobicity[residue]

            if current_score > max_score:
                max_score = current_score
                best_frame = (a_pos_residues, d_pos_residues)
        return best_frame

    def predict_oligomer(a_res, d_res):
        """
        Predicts oligomeric state based on core residue composition rules.
        """
        # Rule 1: Asparagine (N) at a core position is a strong dimer signal.
        if 'N' in a_res or 'N' in d_res:
            return 2

        a_counts = collections.Counter(a_res)
        d_counts = collections.Counter(d_res)

        # Rule 2: Alanine (A) at 'a' and Threonine (T) at 'd' signals a tetramer.
        # This handles imperfect repeats by checking if the majority of residues match.
        if a_counts.get('A', 0) >= len(a_res) * 0.75 and d_counts.get('T', 0) >= len(d_res) * 0.75:
            return 4
        
        # Rule 3: Glutamine (Q) in the core is often found in trimers.
        if 'Q' in a_res or 'Q' in d_res:
            return 3

        # Rule 4: Alanine at 'a' with other large hydrophobes (I, L, W) at 'd' is trimeric.
        # Also Alanine at 'a' and Isoleucine at 'd' can be a trimer.
        if a_counts.get('A', 0) >= len(a_res) * 0.75:
            if d_counts.get('L',0) > 0 or d_counts.get('W',0) > 0 or d_counts.get('I',0) > 0:
                 return 3

        # Default prediction: Trimer is a common state for coiled-coils.
        return 3

    results = []
    for seq in sequences:
        a_residues, d_residues = find_best_frame(seq)
        state = predict_oligomer(a_residues, d_residues)
        results.append(state)
        
    # The problem asks to output the final "equation" showing each number.
    # We will format the output string as requested.
    result_str = ", ".join(map(str, results))
    print(f"The predicted oligomeric states for the sequences are:")
    print(result_str)

solve()