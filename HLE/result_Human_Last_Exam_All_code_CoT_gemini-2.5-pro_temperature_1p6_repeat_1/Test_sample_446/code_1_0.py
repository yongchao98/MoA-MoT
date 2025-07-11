def find_integrin_binding_peptide():
    """
    Identifies the most likely integrin-binding peptide from a list of candidates
    by checking against known high-affinity sequences.
    """

    # List of well-known integrin-binding motifs. RGDSP from fibronectin is a classic
    # example used extensively in in vitro assays.
    known_motifs = ["RGDSP", "RGDV", "LDV"]

    # The candidate peptides from the multiple-choice question.
    candidates = {
        "A": "RGDMAA",
        "B": "RGDSPSS",
        "C": "RGDLTTP",
        "D": "RGDQVSK",
        "E": "RGDARGG"
    }

    print("Searching for known integrin-binding motifs in the candidate peptides...")
    print(f"Known motifs database: {known_motifs}\n")

    found_peptide = None
    reason = ""

    # Iterate through each candidate peptide to see if it contains a known motif.
    # We prioritize the RGDSP motif as it is exceptionally common in research.
    for key, peptide in candidates.items():
        # Check if the peptide starts with the classic "RGDSP" sequence.
        if peptide.startswith("RGDSP"):
            found_peptide = key
            reason = f"Peptide '{peptide}' contains the classic 'RGDSP' motif from fibronectin, which is known to bind integrin receptors with high affinity and is widely used in in vitro assays."
            break # Found the most likely candidate

    # As a secondary check, although RGDLTTP is also a known binder from protein CYR61,
    # RGDSP is arguably more foundational and common in standard assays.
    if not found_peptide:
         for key, peptide in candidates.items():
             if peptide.startswith("RGDLTT"): # Check for the other known binder
                found_peptide = key
                reason = f"Peptide '{peptide}' contains the 'RGDLTT' motif from the protein CYR61, which is also a known integrin ligand."
                break

    if found_peptide:
        print(f"Result: Candidate {found_peptide} ('{candidates[found_peptide]}') is the most likely to bind an integrin receptor.")
        print(f"Reason: {reason}")
    else:
        print("No definitive classic motif found in the primary search of candidates.")

find_integrin_binding_peptide()