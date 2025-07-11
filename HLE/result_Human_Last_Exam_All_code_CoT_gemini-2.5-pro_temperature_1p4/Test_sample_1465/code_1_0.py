import collections

def analyze_coiled_coils(sequences):
    """
    Analyzes a list of coiled-coil sequences to predict their oligomeric state.

    Args:
        sequences (list): A list of protein sequence strings.

    Returns:
        list: A list of predicted integer oligomeric states.
    """
    final_predictions = []
    print("Analyzing Coiled-Coil Sequences...\n")

    for i, seq in enumerate(sequences):
        a_pos = []
        d_pos = []

        # The heptad repeat is 'abcdefg'. We assume the provided sequences start at 'g'.
        # Therefore, 'a' is at index 1, 8, 15... (python index)
        # and 'd' is at index 4, 11, 18... (python index)
        for j in range(len(seq)):
            # The register is offset by 1 because we assume the sequence starts at 'g' (pos 6 of 0-6)
            # a = pos 0, b = pos 1, ..., g = pos 6
            # So, (j+1) % 7 gives the position in the 'gabcdef' frame.
            # We look for 'a' (pos 1) and 'd' (pos 4) in the 'abcdefg' frame
            pos_in_heptad = j % 7
            if pos_in_heptad == 1:  # 'a' position
                a_pos.append(seq[j])
            elif pos_in_heptad == 4:  # 'd' position
                d_pos.append(seq[j])
        
        a_counts = collections.Counter(a_pos)
        d_counts = collections.Counter(d_pos)

        prediction = 0
        reason = ""

        # Rule 1: Tetramer (4-mer) if 'a' and 'd' are dominated by Isoleucine (I)
        if a_counts['I'] >= len(a_pos) - 1 and d_counts['I'] >= len(d_pos) - 2:
            prediction = 4
            reason = "Core dominated by Isoleucine (I) at 'a' and 'd' positions favors a tetramer."
        
        # Rule 2: Trimer (3-mer) if 'a' is Isoleucine (I) and 'd' is Threonine (T)
        elif a_counts['I'] >= len(a_pos) -1 and d_counts['T'] >= len(d_pos) - 2:
            prediction = 3
            reason = "Core with Isoleucine (I) at 'a' and Threonine (T) at 'd' is a classic trimer motif."

        # Rule 3: Dimer (2-mer) if a polar Asparagine (N) is in the core 'a' position
        elif 'N' in a_counts:
            prediction = 2
            reason = "Polar Asparagine (N) at an 'a' position is a hallmark of a parallel dimer."
        
        # Rule 4: Default to Dimer (2-mer) for other cases
        else:
            prediction = 2
            reason = "No strong trimer/tetramer motif detected. Defaulting to dimer, the most common state."

        final_predictions.append(prediction)
        
        print(f"Sequence {i+1}: {seq}")
        print(f" -> 'a' positions: {a_pos}")
        print(f" -> 'd' positions: {d_pos}")
        print(f" -> Prediction: {prediction}-mer")
        print(f" -> Reason: {reason}\n")

    print("="*30)
    print("Final predicted oligomeric states:")
    # The required final output format is each number in the equation.
    print(*final_predictions, sep=",")
    print("="*30)


# The protein sequences provided by the user
protein_sequences = [
    "EIAQALKEIAKALKEIAWALKEIAQALK",
    "EIAALKQEIAALKKENAALKQEIAALKQ",
    "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
    "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
    "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
]

analyze_coiled_coils(protein_sequences)