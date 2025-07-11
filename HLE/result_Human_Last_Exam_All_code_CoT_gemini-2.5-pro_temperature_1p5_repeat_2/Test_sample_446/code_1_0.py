def find_integrin_binding_peptide():
    """
    Identifies the most likely integrin-binding peptide from a list of candidates.
    """
    # Dictionary of peptide choices
    peptides = {
        "A": "RGDMAA",
        "B": "RGDSPSS",
        "C": "RGDLTTP",
        "D": "RGDQVSK",
        "E": "RGDARGG"
    }

    # A mock database of a well-known, high-affinity integrin-binding motif.
    # The RGDSP sequence is from fibronectin and is a canonical integrin ligand.
    known_motif = "RGDSP"

    print(f"Searching for the canonical integrin-binding motif '{known_motif}' in the following peptides:")
    for key, peptide in peptides.items():
        print(f"Choice {key}: {peptide}")

    print("\n--- Analysis ---")
    found_peptide = None
    found_key = None

    for key, peptide in peptides.items():
        if known_motif in peptide:
            found_peptide = peptide
            found_key = key
            break

    if found_peptide:
        print(f"Result: Peptide '{found_peptide}' (Choice {found_key}) contains the known high-affinity motif '{known_motif}'.")
        print("Explanation: While all peptides have the RGD core, the 'SP' (Serine-Proline) dipeptide following RGD is well-documented to induce a beta-turn conformation.")
        print("This specific structure presents the RGD sequence optimally for binding to integrin receptors, making it a very common peptide used in in vitro assays.")
    else:
        print("No peptide containing the canonical RGDSP motif was found.")

find_integrin_binding_peptide()