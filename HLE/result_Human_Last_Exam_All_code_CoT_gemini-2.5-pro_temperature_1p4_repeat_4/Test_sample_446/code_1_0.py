def find_integrin_binder():
    """
    Identifies the most likely integrin-binding peptide from a list of choices.

    Integrin receptors often bind to the RGD motif. The flanking sequences are
    critical for binding affinity. The sequence RGDSP is a well-characterized
    motif from the natural ligand fibronectin, known for its high-affinity
    binding to integrins like alpha-v beta-3 and alpha-5 beta-1. This script
    identifies the peptide containing this canonical sequence.
    """
    peptides = {
        "A": "RGDMAA",
        "B": "RGDSPSS",
        "C": "RGDLTTP",
        "D": "RGDQVSK",
        "E": "RGDARGG"
    }

    # The RGDSP motif is a known high-affinity binder from fibronectin.
    # We will search for the peptide containing this motif.
    high_affinity_motif = "RGDSP"
    best_candidate = None
    best_candidate_choice = None

    for choice, peptide in peptides.items():
        if high_affinity_motif in peptide:
            best_candidate = peptide
            best_candidate_choice = choice
            break

    if best_candidate:
        print(f"Analysis of Peptides:")
        for choice, peptide in peptides.items():
            print(f"  {choice}. {peptide}")
        
        print(f"\nConclusion:")
        print(f"The peptide motif '{high_affinity_motif}' is a well-known sequence from fibronectin, a natural ligand for integrins.")
        print("This motif confers high-affinity binding.")
        print(f"The peptide containing this sequence is '{best_candidate}'.")
        print(f"Therefore, choice {best_candidate_choice} is the most likely to bind an integrin receptor.")

find_integrin_binder()