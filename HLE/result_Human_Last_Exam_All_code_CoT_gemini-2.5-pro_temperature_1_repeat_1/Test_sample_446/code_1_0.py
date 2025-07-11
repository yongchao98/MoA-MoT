def find_integrin_binding_peptide():
    """
    Identifies a likely integrin-binding peptide from a list of candidates.

    The RGD (Arginine-Glycine-Aspartic acid) sequence is a primary recognition
    site for many integrin receptors. However, the flanking amino acids are crucial
    for specificity and high-affinity binding.

    One of the most well-studied integrin-binding motifs is RGDSP, derived from
    fibronectin, a major extracellular matrix protein. This script checks which of
    the candidate peptides contains this well-known motif.
    """
    
    # List of candidate peptides from the problem
    candidate_peptides = {
        "A": "RGDMAA",
        "B": "RGDSPSS",
        "C": "RGDLTTP",
        "D": "RGDQVSK",
        "E": "RGDARGG"
    }

    # A small, simulated database of well-known integrin-binding sequences
    known_binders = ["RGDSP"]

    print("Searching for the most likely integrin-binding peptide...")
    print("Candidate Peptides:")
    for key, peptide in candidate_peptides.items():
        print(f"{key}: {peptide}")
    
    print("\nKnown high-affinity motif being searched for:", known_binders[0])
    
    found_peptide = None
    found_key = None

    # Search for a candidate that contains a known binding motif
    for key, peptide in candidate_peptides.items():
        for binder in known_binders:
            if binder in peptide:
                found_peptide = peptide
                found_key = key
                break
        if found_peptide:
            break
            
    if found_peptide:
        print(f"\nResult: Peptide '{found_peptide}' (Choice {found_key}) contains the known integrin-binding motif '{known_binders[0]}'.")
        print("This sequence is derived from fibronectin and is well-documented to bind integrin receptors in vitro.")
    else:
        print("\nResult: No peptide containing the specified known motifs was found.")

find_integrin_binding_peptide()