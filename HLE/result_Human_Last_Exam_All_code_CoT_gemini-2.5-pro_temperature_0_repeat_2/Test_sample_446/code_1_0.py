def find_integrin_binding_peptide():
    """
    Identifies the most likely integrin-binding peptide from a list of choices
    based on known binding motifs.
    """
    # List of candidate peptides from the answer choices
    peptides = {
        "A": "RGDMAA",
        "B": "RGDSPSS",
        "C": "RGDLTTP",
        "D": "RGDQVSK",
        "E": "RGDARGG"
    }

    # A small database of well-known integrin-binding motifs.
    # The RGDSP sequence from fibronectin is a canonical example.
    known_motifs = ["RGDSP"]

    print("Analyzing peptides for known integrin-binding motifs...")
    print("-" * 30)

    best_candidate = None
    best_candidate_key = None

    for key, peptide in peptides.items():
        print(f"Checking Candidate {key}: {peptide}")
        for motif in known_motifs:
            if peptide.startswith(motif):
                best_candidate = peptide
                best_candidate_key = key
                print(f"-> Found a match! The peptide contains the known motif '{motif}'.")
                break
        if best_candidate:
            break

    print("-" * 30)
    if best_candidate:
        print("Conclusion:")
        print("The RGD (Arginine-Glycine-Aspartic acid) sequence is the primary recognition site for many integrins.")
        print("However, the flanking amino acids are critical for binding affinity and specificity.")
        print(f"The sequence '{known_motifs[0]}' is a well-characterized motif from fibronectin, known to bind potently to integrins like α5β1.")
        print(f"Peptide {best_candidate_key}, '{best_candidate}', contains this canonical motif and is therefore the most likely to bind an integrin receptor.")
    else:
        print("No peptide containing a known high-affinity motif was found in the database.")

    # The final answer is the letter corresponding to the best candidate.
    if best_candidate_key:
        print(f"\nFinal Answer: {best_candidate_key}")


find_integrin_binding_peptide()