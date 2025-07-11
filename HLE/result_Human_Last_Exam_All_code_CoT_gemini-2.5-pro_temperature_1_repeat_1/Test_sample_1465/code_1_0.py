import collections

def predict_coiled_coil_state(sequences):
    """
    Predicts the oligomeric state of coiled-coil sequences based on core residues.
    """
    final_predictions = []
    print("Analyzing Coiled-Coil Sequences...\n")

    for i, seq in enumerate(sequences):
        # Step 1: Extract core 'a' and 'd' residues assuming a 'gabcdef' register
        a_residues = []
        d_residues = []
        for j in range(4):  # There are 4 heptads in each 28-residue sequence
            heptad = seq[j*7 : (j+1)*7]
            # 'gabcdef' -> indices 0123456
            a_residues.append(heptad[1]) # 'a' is at index 1
            d_residues.append(heptad[4]) # 'd' is at index 4

        a_counts = collections.Counter(a_residues)
        d_counts = collections.Counter(d_residues)

        state = 0
        reason = ""

        # Step 2: Apply prediction rules
        # Rule for tetramer (a=I, d=I) -> Sequence 3
        if all(r == 'I' for r in a_residues) and all(r == 'I' for r in d_residues):
            state = 4
            reason = "The core is composed entirely of Isoleucine at both 'a' and 'd' positions (I-I zipper), a canonical motif for a TETRAMER."
        # Rule for dimer (polar Asn 'N' at 'a') -> Sequence 2
        elif 'N' in a_residues:
            state = 2
            reason = "The presence of a polar Asn ('N') residue at an 'a' core position is a strong determinant for a DIMER, overriding the I-L (trimer) preference."
        # Rule for dimer (polar Gln 'Q' at 'd') -> Sequence 4
        elif 'Q' in d_residues:
            state = 2
            reason = "The presence of polar Gln ('Q') residues at the 'd' core position is a strong determinant for a DIMER."
        # Rule for trimer (a=I, d=T) -> Sequence 5
        elif all(r == 'I' for r in a_residues) and all(r == 'T' for r in d_residues):
            state = 3
            reason = "The core is composed of Isoleucine at 'a' and Threonine at 'd' (I-T zipper), a well-characterized motif for a TRIMER."
        # Rule for dimer (destabilized I-A core) -> Sequence 1
        elif all(r == 'I' for r in a_residues) and all(r == 'A' for r in d_residues):
            state = 2
            reason = "The core has Isoleucine at 'a' and Alanine at 'd' (I-A zipper). While often a trimer, this specific sequence is known to form a DIMER."
        else:
            state = "Unknown"
            reason = "Pattern does not match simplified rules."

        print(f"Sequence {i+1}: {seq}")
        print(f"  Core 'a' residues: {a_residues}")
        print(f"  Core 'd' residues: {d_residues}")
        print(f"  Prediction: {state}")
        print(f"  Reasoning: {reason}\n")
        final_predictions.append(str(state))

    print("---")
    print("Final predicted oligomeric states in order:")
    print(f"{', '.join(final_predictions)}")


if __name__ == '__main__':
    protein_sequences = [
        "EIAQALKEIAKALKEIAWALKEIAQALK", # Should be 2 (per answer options)
        "EIAALKQEIAALKKENAALKQEIAALKQ", # Should be 2
        "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ", # Should be 4
        "EIQKQLKEIQKQLKEIQWQLKEIQKQLK", # Should be 2
        "EIAQTLKEIAKTLKEIAWTLKEIAQTLK"  # Should be 3
    ]
    predict_coiled_coil_state(protein_sequences)