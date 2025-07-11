def find_integrin_binding_peptide():
    """
    Analyzes a list of RGD-containing peptides to find the one most likely
    to bind an integrin receptor based on known high-affinity motifs.
    """
    # A dictionary of the candidate peptides provided in the question
    peptides = {
        "A": "RGDMAA",
        "B": "RGDSPSS",
        "C": "RGDLTTP",
        "D": "RGDQVSK",
        "E": "RGDARGG"
    }

    # A mock database of well-known high-affinity integrin-binding motifs.
    # The RGDSP sequence is a famous motif from fibronectin, a natural ligand.
    known_high_affinity_motifs = {
        "RGDSP": "a classic high-affinity binding site from the protein fibronectin, which binds strongly to integrins like α5β1 and αvβ3."
    }

    print("Analyzing candidate peptides against a database of known motifs...")
    print("-" * 60)

    found_candidate = False
    for option, peptide_sequence in peptides.items():
        for motif, description in known_high_affinity_motifs.items():
            # Check if the known high-affinity motif is present in the candidate peptide
            if motif in peptide_sequence:
                print(f"Result: Option {option} is the most likely candidate.")
                print(f"Peptide Sequence: {peptide_sequence}")
                print(f"Reason: This peptide contains the '{motif}' sequence. This is {description}")
                found_candidate = True
                break  # Exit the inner loop once a match is found
        if found_candidate:
            break # Exit the outer loop as we've found our answer

    if not found_candidate:
        print("None of the peptides contain a known high-affinity motif from the database.")
        print("While all contain the core RGD sequence, further analysis would be required.")

    print("-" * 60)

# Run the analysis function
find_integrin_binding_peptide()