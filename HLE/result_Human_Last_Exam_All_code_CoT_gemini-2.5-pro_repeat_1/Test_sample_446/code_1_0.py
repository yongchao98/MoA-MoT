def find_integrin_binding_peptide():
    """
    Analyzes a list of RGD-containing peptides to find the one most likely
    to bind an integrin receptor based on known motifs.

    The RGDSP sequence, found in fibronectin, is a well-established
    high-affinity binding motif for certain integrins (e.g., a5b1). This script
    checks which of the candidate peptides contains this known motif.
    """
    
    # The list of peptide options
    peptides = {
        "A": "RGDMAA",
        "B": "RGDSPSS",
        "C": "RGDLTTP",
        "D": "RGDQVSK",
        "E": "RGDARGG"
    }

    # A known high-affinity integrin binding motif derived from fibronectin
    known_motif = "RGDSP"
    
    most_likely_choice = None
    
    print(f"Analyzing peptides for the known high-affinity motif: '{known_motif}'")
    print("---------------------------------------------------------")

    for choice, sequence in peptides.items():
        # Check if the peptide's sequence starts with the known motif
        if sequence.startswith(known_motif):
            print(f"Peptide {choice} ('{sequence}') contains the '{known_motif}' motif. It is highly likely to bind an integrin.")
            most_likely_choice = choice
        else:
            print(f"Peptide {choice} ('{sequence}') does not contain the specific '{known_motif}' motif.")

    print("---------------------------------------------------------")
    if most_likely_choice:
        print(f"Conclusion: The peptide most likely to bind an integrin receptor is from choice {most_likely_choice}.")
    else:
        print("Conclusion: No peptide matching the specific known motif was found.")

find_integrin_binding_peptide()