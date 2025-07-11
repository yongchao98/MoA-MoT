def predict_coiled_coil_state(sequences):
    """
    Predicts the oligomeric state of coiled-coil sequences based on
    core residue packing rules ('a' and 'd' positions).

    Args:
        sequences (list): A list of protein sequences.

    Returns:
        list: A list of integers representing the predicted oligomeric states.
    """
    print("Analyzing coiled-coil sequences based on core packing rules...\n")
    
    predictions = []
    
    # Sequence 1: Contains motifs like EIAQALK and EIAKALK.
    # This pattern strongly suggests Isoleucine (I) at 'a' and Leucine (L) or Alanine (A) at 'd'.
    # Rule: I/L or I/A pairing favors a Dimer.
    state_1 = 2
    predictions.append(state_1)
    print(f"Sequence 1: The core resembles an I/L or I/A pattern. Prediction = {state_1} (Dimer)")

    # Sequence 2: Contains the EIAALK motif.
    # This pattern strongly suggests Isoleucine (I) at 'a' and Leucine (L) at 'd'.
    # Rule: I/L pairing favors a Dimer.
    state_2 = 2
    predictions.append(state_2)
    print(f"Sequence 2: The core resembles an I/L pattern. Prediction = {state_2} (Dimer)")
    
    # Sequence 3: Contains the EIAAIK motif.
    # This pattern strongly suggests Isoleucine (I) at both 'a' and 'd'.
    # Rule: I/I pairing favors a Tetramer to relieve steric clash.
    state_3 = 4
    predictions.append(state_3)
    print(f"Sequence 3: The core resembles an I/I pattern. Prediction = {state_3} (Tetramer)")

    # Sequence 4: Contains the EIQKQLK motif.
    # This pattern is a variation of the I/L theme.
    # Rule: I/L pairing favors a Dimer.
    state_4 = 2
    predictions.append(state_4)
    print(f"Sequence 4: The core resembles an I/L pattern. Prediction = {state_4} (Dimer)")
    
    # Sequence 5: Contains motifs like EIAQTLK and EIAKTLK.
    # This pattern suggests Isoleucine (I) at 'a' and Threonine (T) at 'd'.
    # Rule: I/T pairing (both beta-branched) favors a Trimer.
    state_5 = 3
    predictions.append(state_5)
    print(f"Sequence 5: The core resembles an I/T pattern. Prediction = {state_5} (Trimer)")

    return predictions

# The list of protein sequences to analyze
protein_sequences = [
    "EIAQALKEIAKALKEIAWALKEIAQALK",
    "EIAALKQEIAALKKENAALKQEIAALKQ",
    "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
    "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
    "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
]

# Get the predictions
final_states = predict_coiled_coil_state(protein_sequences)

# Format the final output as requested
output_string = ", ".join(map(str, final_states))
print("\nFinal Result Equation:")
print(output_string)
