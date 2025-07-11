import sys

def find_best_integrin_binder():
    """
    Analyzes a list of RGD-containing peptides to identify the one most likely to bind
    to an integrin receptor based on known high-affinity motifs.
    """
    # List of peptide candidates from the answer choices
    peptides = {
        'A': 'RGDMAA',
        'B': 'RGDSPSS',
        'C': 'RGDLTTP',
        'D': 'RGDQVSK',
        'E': 'RGDARGG'
    }

    # Known high-affinity motif. The RGDSP sequence is from fibronectin.
    # The Ser-Pro (SP) pair is known to promote a beta-turn structure,
    # which is critical for high-affinity binding to many integrins.
    known_motif = "RGDSP"
    best_candidate_letter = None
    best_candidate_peptide = None

    print("Analyzing peptides for integrin binding potential...\n")

    for letter, peptide in peptides.items():
        # Check if the peptide contains the known high-affinity motif
        if known_motif in peptide:
            print(f"Candidate {letter}: {peptide} - Contains the known high-affinity '{known_motif}' motif from fibronectin.")
            print("   -> This sequence is highly likely to form a stable structure for integrin binding.")
            best_candidate_letter = letter
            best_candidate_peptide = peptide
        else:
            print(f"Candidate {letter}: {peptide} - Does not contain the classic high-affinity '{known_motif}' motif.")
    
    if best_candidate_letter:
        print("\n--- Conclusion ---")
        print(f"The peptide '{best_candidate_peptide}' (Choice {best_candidate_letter}) is the most likely candidate to bind an integrin receptor in an in vitro assay.")
        print("This is because it contains the RGDSP sequence, a well-characterized motif from fibronectin known to promote a specific conformation that leads to high-affinity binding.")
        # To satisfy the prompt to output the final equation/answer format
        sys.stdout.write(f'<<<B>>>')
    else:
        print("\nNo peptide with the primary known high-affinity motif was found.")
        sys.stdout.write(f'<<<No clear answer found>>>')

# Execute the analysis
find_best_integrin_binder()