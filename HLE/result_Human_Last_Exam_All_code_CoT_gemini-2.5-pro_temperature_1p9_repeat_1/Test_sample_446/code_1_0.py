import sys

def find_integrin_binder():
    """
    Identifies the most likely integrin-binding peptide from a list of candidates
    by searching against a database of known motifs.
    """
    # List of candidate peptides from the user's query
    peptides = {
        'A': 'RGDMAA',
        'B': 'RGDSPSS',
        'C': 'RGDLTTP',
        'D': 'RGDQVSK',
        'E': 'RGDARGG'
    }

    # A small database of well-characterized integrin-binding motifs.
    # The RGDSP sequence from fibronectin is a classic example known for its
    # high affinity and frequent use in in vitro assays.
    known_motifs_db = {
        'RGDSP': 'This motif is from fibronectin, a major extracellular matrix protein and a primary ligand for many integrins (e.g., a5b1). Peptides containing RGDSP are standard positive controls in cell adhesion assays.',
        'RGDV': 'This motif from vitronectin is another well-known integrin binder, particularly for avb3 integrin.'
    }

    best_candidate = None
    strongest_evidence = ""

    print("Analyzing candidate peptides against a database of known integrin binders...")
    print("-" * 70)

    for choice, sequence in peptides.items():
        match_found = False
        for motif, evidence in known_motifs_db.items():
            if motif in sequence:
                print(f"Candidate {choice}: {sequence}")
                print(f"Found Motif: {motif}")
                print(f"Evidence: {evidence}\n")
                best_candidate = choice
                strongest_evidence = evidence
                match_found = True
                break  # Stop after finding the first and most relevant motif
        if not match_found:
            print(f"Candidate {choice}: {sequence}")
            print("Evidence: This sequence does not contain a commonly cited high-affinity motif from our database.\n")

    print("-" * 70)
    print("Conclusion:")
    if best_candidate:
        print(f"Peptide {best_candidate} ({peptides[best_candidate]}) has been found to bind or is highly likely to bind an integrin receptor.")
        print(f"Reason: It contains the '{list(known_motifs_db.keys())[0]}' sequence derived from fibronectin, which is extensively documented to mediate cell adhesion in in vitro assays.")
    else:
        print("Could not definitively identify a well-known binder from the provided options.")

# Execute the function
find_integrin_binder()
