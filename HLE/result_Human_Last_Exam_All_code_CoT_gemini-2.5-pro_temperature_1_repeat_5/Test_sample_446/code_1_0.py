def find_most_likely_binder():
    """
    Identifies the peptide most likely to bind to an integrin receptor
    by comparing it to a known high-affinity sequence from fibronectin.
    """
    # All choices contain the core RGD motif. Binding affinity is modulated
    # by the flanking amino acids.
    candidates = {
        "A": "RGDMAA",
        "B": "RGDSPSS",
        "C": "RGDLTTP",
        "D": "RGDQVSK",
        "E": "RGDARGG"
    }

    # The sequence from fibronectin, a well-known integrin ligand, is RGDSPASS.
    # Peptides mimicking this sequence show high affinity.
    reference_sequence = "RGDSPASS"

    best_match_label = None
    highest_score = -1

    print("Analyzing peptide candidates based on similarity to the fibronectin binding site (RGDSPASS):\n")

    # Calculate a similarity score for each candidate
    for label, peptide in candidates.items():
        score = 0
        # Compare each amino acid position-wise
        for i in range(min(len(peptide), len(reference_sequence))):
            if peptide[i] == reference_sequence[i]:
                score += 1
        
        print(f"Candidate {label} ({peptide}): Similarity Score = {score}")

        if score > highest_score:
            highest_score = score
            best_match_label = label

    best_peptide = candidates[best_match_label]
    
    print("\n---")
    print(f"Conclusion: The peptide '{best_peptide}' (Choice {best_match_label}) has the highest similarity score.")
    print("This sequence is very similar to the natural fibronectin sequence RGDSPASS, making it the most likely candidate to bind an integrin receptor in an in vitro assay.")

find_most_likely_binder()