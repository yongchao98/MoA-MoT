def find_integrin_binding_peptide():
    """
    Identifies the most likely integrin-binding peptide from a list of candidates.

    This function is based on the biological knowledge that while the RGD motif is essential,
    flanking sequences greatly influence binding affinity. The RGDSP sequence is a
    well-characterized, high-affinity motif found in fibronectin, a major ligand for integrins.
    """

    # List of candidate peptides from the answer choices
    candidates = {
        "A": "RGDMAA",
        "B": "RGDSPSS",
        "C": "RGDLTTP",
        "D": "RGDQVSK",
        "E": "RGDARGG"
    }

    # A knowledge base of well-established high-affinity integrin-binding motifs.
    # The RGDSP motif is a classic example.
    known_high_affinity_motifs = ["RGDSP"]

    # The best candidate found so far
    best_candidate_key = None
    best_candidate_peptide = None
    reason = "No peptide contains a known high-affinity motif."

    # Iterate through each candidate peptide
    for key, peptide in candidates.items():
        # Check if the peptide contains a known high-affinity motif
        for motif in known_high_affinity_motifs:
            if motif in peptide:
                best_candidate_key = key
                best_candidate_peptide = peptide
                reason = (f"Peptide {peptide} contains the well-known high-affinity "
                          f"integrin-binding motif '{motif}', which is derived from fibronectin.")
                break # Found a match, no need to check other motifs for this peptide
        if best_candidate_key:
            break # Found the best candidate, stop searching

    # Print the result
    if best_candidate_key:
        print(f"Answer Choice: {best_candidate_key}")
        print(f"Peptide: {best_candidate_peptide}")
        print(f"Reasoning: {reason}")
    else:
        # This part will not be reached for the given problem, but is good practice.
        print("Could not identify a clear winner based on the known motifs.")
        print(f"All peptides contain the base 'RGD' motif, but an additional well-known "
              f"high-affinity sequence like 'RGDSP' was not found in any candidate.")


if __name__ == "__main__":
    find_integrin_binding_peptide()