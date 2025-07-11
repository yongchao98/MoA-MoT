import collections

def predict_oligomeric_state(sequence):
    """
    Predicts the oligomeric state of a coiled-coil sequence based on a simplified
    analysis of residues at the 'a' and 'd' positions of the heptad repeat.
    """
    a_positions = []
    d_positions = []

    # Assume standard heptad repeat starting from the first residue
    # and extract amino acids at positions 'a' (index 0) and 'd' (index 3)
    for i in range(0, len(sequence), 7):
        heptad = sequence[i:i+7]
        if len(heptad) >= 4:
            a_positions.append(heptad[0])
            d_positions.append(heptad[3])

    a_counts = collections.Counter(a_positions)
    d_counts = collections.Counter(d_positions)
    total_heptads = len(a_positions)

    # Rule 1: N at a core 'a' position strongly promotes trimers.
    if 'N' in a_positions:
        return 3, "Contains 'N' at a core 'a' position, which favors trimers."

    # Rule 2: 'I' at 'a' and 'A' at 'd' is characteristic of tetramers.
    if a_counts.get('I', 0) / total_heptads > 0.8 and d_counts.get('A', 0) / total_heptads > 0.8:
        return 4, "Core with 'I' at 'a' and 'A' at 'd' favors tetramers."

    # Rule 3: 'I' at 'a' and 'L' at 'd' is the classic dimer signature.
    if a_counts.get('I', 0) / total_heptads > 0.8 and d_counts.get('L', 0) / total_heptads > 0.8:
        return 2, "Core with 'I' at 'a' and 'L' at 'd' strongly favors dimers."

    # Rule 4: A mixed, destabilized core with a bulky 'W' can favor a trimer.
    if 'W' in d_positions and ('K' in d_positions or 'Q' in d_positions):
        return 3, "A mixed core with a bulky 'W' at 'd' can be accommodated by a trimer."

    # Rule 5: A polar core with 'Q' at 'd' can be accommodated in a dimer.
    if d_counts.get('Q', 0) / total_heptads > 0.6:
        return 2, "A polar core with 'Q' at 'd' is compatible with a dimer."

    return "Unknown", "Sequence did not match a clear prediction rule."

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

    print("--- Oligomeric State Prediction ---")
    predictions = []
    for i, seq in enumerate(sequences):
        state, reason = predict_oligomeric_state(seq)
        print(f"\nSequence {i+1}: {seq}")
        print(f"Prediction: {state}-mer")
        print(f"Reason: {reason}")
        predictions.append(state)

    print("\n--- Final Result ---")
    print(f"The predicted oligomeric states for the sequences are: {', '.join(map(str, predictions))}")
    # The final equation as requested.
    print(f"{predictions[0]}, {predictions[1]}, {predictions[2]}, {predictions[3]}, {predictions[4]}.")


if __name__ == "__main__":
    main()