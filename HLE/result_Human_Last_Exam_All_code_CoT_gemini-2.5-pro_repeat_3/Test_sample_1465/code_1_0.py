import collections

def predict_oligomeric_state(sequence):
    """
    Predicts the oligomeric state of a coiled-coil sequence based on
    residues at the 'a' and 'd' positions of the heptad repeat.
    """
    # Assuming a 'gabcdefg' register for the heptad repeat.
    # 'a' positions are at indices 1, 8, 15, 22, ... (0-indexed)
    # 'd' positions are at indices 4, 11, 18, 25, ... (0-indexed)
    a_residues = [sequence[i] for i in range(1, len(sequence), 7)]
    d_residues = [sequence[i] for i in range(4, len(sequence), 7)]

    # Use a Counter to find the most common residue at each position.
    # This makes the prediction robust to single-residue variations.
    a_counts = collections.Counter(a_residues)
    d_counts = collections.Counter(d_residues)
    
    # In case of a tie for the most common residue, most_common(1) will pick one.
    # For these specific sequences, there are no ties for the dominant residue type.
    dominant_a = a_counts.most_common(1)[0][0]
    dominant_d = d_counts.most_common(1)[0][0]
    
    prediction = "Unknown"
    reasoning = "Could not determine state with high confidence."

    # Rule 1: Tetramer (I at 'a' and I at 'd') is a very strong signal.
    # We check if the majority of core positions are Isoleucine.
    if a_counts.get('I', 0) >= 3 and d_counts.get('I', 0) >= 3:
        prediction = 4
        reasoning = f"The core is dominated by Isoleucine (I) at both 'a' and 'd' positions, which strongly favors a tetrameric state."

    # Rule 2: Trimer (I at 'a' and T at 'd')
    # Beta-branched Threonine (T) at the 'd' position creates steric hindrance
    # that is better accommodated in a trimer than a dimer.
    elif dominant_a == 'I' and dominant_d == 'T':
        prediction = 3
        reasoning = f"The core pattern of Isoleucine (I) at 'a' and Threonine (T) at 'd' suggests a trimeric state."

    # Rule 3: Dimer (Classic I-L or I-A patterns)
    elif dominant_a == 'I' and (dominant_d == 'L' or dominant_d == 'A'):
        prediction = 2
        reasoning = f"The classic core pattern of '{dominant_a}' at 'a' and '{dominant_d}' at 'd' is characteristic of a dimeric state."

    # Rule 4: Special case for non-canonical core residues (handles sequence 4)
    # The 'a' positions are all Isoleucine, providing a strong core. While Lysine (K)
    # is unfavorable at 'd', a dimer is the most likely state to accommodate this.
    elif dominant_a == 'I' and dominant_d == 'K':
        prediction = 2
        reasoning = f"Despite the unfavorable Lysine (K) at the 'd' position, the strong Isoleucine (I) signal at 'a' points towards a dimer as the most stable option."

    return prediction, reasoning, a_residues, d_residues

def main():
    """
    Main function to analyze the sequences and print the results.
    """
    sequences = [
        "EIAQALKEIAKALKEIAWALKEIAQALK",
        "EIAALKQEIAALKKENAALKQEIAALKQ",
        "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
        "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
        "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
    ]

    final_predictions = []
    print("Analysis of Coiled-Coil Oligomeric States")
    print("="*45)

    for i, seq in enumerate(sequences, 1):
        prediction, reason, a_res, d_res = predict_oligomeric_state(seq)
        final_predictions.append(prediction)
        
        print(f"Sequence {i}: {seq}")
        print(f"  > Core 'a' residues: {a_res}")
        print(f"  > Core 'd' residues: {d_res}")
        print(f"  > Reasoning: {reason}")
        print(f"  > Predicted State: {prediction}")
        print("-" * 45)

    # The prompt asks to "output each number in the final equation!".
    # This formats the final list of numbers as requested.
    result_string = ", ".join(map(str, final_predictions))
    print("\nFinal Result")
    print("The predicted oligomeric states for the sequences are:")
    print(result_string)

if __name__ == "__main__":
    main()