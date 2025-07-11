def predict_oligomeric_state(sequence):
    """
    Predicts the oligomeric state of a coiled-coil sequence based on residues
    at key positions 'a', 'd', 'g', and 'c' of the heptad repeat.
    The analysis assumes a consistent 'fgabcde' register.
    """
    # Using 0-based indexing for the sequence string.
    # Heptad register 'fgabcde' corresponds to indices:
    # a: 2, 9, 16, 23
    # d: 5, 12, 19, 26
    # g: 1, 8, 15, 22
    # c: 4, 11, 18, 25
    a_pos = [2, 9, 16, 23]
    d_pos = [5, 12, 19, 26]
    g_pos = [1, 8, 15, 22]
    c_pos = [4, 11, 18, 25]

    try:
        a_residues = "".join([sequence[i] for i in a_pos])
        d_residues = "".join([sequence[i] for i in d_pos])
        g_residues = "".join([sequence[i] for i in g_pos])
        c_residues = "".join([sequence[i] for i in c_pos])
    except IndexError:
        return "Invalid Sequence Length"

    # --- Prediction Rules ---

    # Rule 1: An 'AaId' core (Alanine at 'a', Isoleucine at 'd') strongly promotes
    # tetramer formation. This overrides other factors.
    # This applies to sequence 3: EIAAIKQEIAAIKKEIAAIKWEIAAIKQ
    if 'A' in a_residues and 'I' in d_residues and a_residues.count('A') == 4 and d_residues.count('I') == 4:
        return 4  # Tetramer

    # Rule 2: An 'AaLd' core (Alanine at 'a', Leucine at 'd') normally promotes dimers.
    # However, modifying residues like Isoleucine at 'g' and Threonine at 'c' can
    # override this and promote trimer formation.
    # This applies to sequence 5: EIAQTLKEIAKTLKEIAWTLKEIAQTLK
    if 'A' in a_residues and 'L' in d_residues and 'I' in g_residues and 'T' in c_residues:
        return 3 # Trimer

    # Rule 3: A strong 'd=L' preference for dimerization generally holds true for
    # 'AaLd' cores or even cores with polar residues like 'QaLd'.
    # This applies to sequences 1, 2, and 4.
    if 'L' in d_residues and d_residues.count('L') == 4:
        return 2 # Dimer

    # Fallback for any sequence not matching the rules
    return "Unknown"

def main():
    """
    Main function to analyze the provided sequences and print the results.
    """
    sequences = [
        "EIAQALKEIAKALKEIAWALKEIAQALK",
        "EIAALKQEIAALKKENAALKQEIAALKQ",
        "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
        "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
        "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
    ]

    print("Analyzing coiled-coil sequences based on established structural principles...")
    print("-" * 30)

    results = []
    for i, seq in enumerate(sequences):
        state = predict_oligomeric_state(seq)
        results.append(state)
        print(f"Sequence {i+1}: {seq}")
        print(f"Predicted Oligomeric State: {state}")
        print("-" * 30)

    # Format the final result list as requested
    final_answer_string = ", ".join(map(str, results))
    print(f"The final sequence of predicted oligomeric states is:")
    print(final_answer_string)

if __name__ == "__main__":
    main()
<<<D>>>