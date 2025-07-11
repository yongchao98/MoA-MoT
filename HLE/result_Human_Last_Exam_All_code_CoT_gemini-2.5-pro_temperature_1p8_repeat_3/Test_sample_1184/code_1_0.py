def find_calcium_channel_hotspots():
    """
    Identifies and prints key residues in the human voltage-gated calcium
    channel beta-1 subunit (CACNB1) based on established structural and
    functional data.

    The residue numbering corresponds to the canonical human CACNB1 sequence
    (UniProt Accession: P54283).
    """

    # Data based on structural analysis (e.g., PDB: 5V1P) and mutagenesis studies.
    # The primary interaction site is the Beta Interaction Domain (BID) on the
    # beta-1 subunit, which forms a deep pocket to accept the Alpha Interaction
    # Domain (AID) from the alpha-1 subunit.

    interaction_hotspots = {
        'Y372': 'Forms a key part of the hydrophobic binding pocket for the alpha-1 AID. Its aromatic ring interacts directly with the AID peptide.',
        'W391': 'A critical residue at the bottom of the binding pocket. The indole ring of Tryptophan 391 makes extensive hydrophobic contact with the AID.',
        'E413': 'Forms an important salt bridge with a conserved arginine residue on the alpha-1 AID, anchoring it in the binding site.',
        'I415': 'Contributes to the hydrophobic floor of the binding pocket, stabilizing the interaction with the AID.'
    }

    # Gating modulation is a direct consequence of the binding event.
    # Therefore, the interaction hotspots are the primary modulatory residues.
    # Additional residues in other domains can fine-tune gating properties.

    gating_modulation_hotspots = {
        'Y372, W391, E413, I415': 'These primary binding residues in the GK domain are the most critical determinants of gating modulation. The binding event itself causes the allosteric changes that affect channel opening, closing, and inactivation.',
        'W300': 'Located in the SH3 domain. This domain can engage in secondary interactions that fine-tune gating, particularly inactivation kinetics.',
        'G302': 'Also in the SH3 domain ligand-binding surface. Mutations in this area can alter the modulatory effect of the beta-1 subunit without completely disrupting the primary interaction.'
    }

    print("Analyzing Human Calcium Channel Beta-1 Subunit (CACNB1, UniProt: P54283)...")
    print("-" * 70)

    # 1. Answer for interaction hotspots
    print("\n[Question 1] Hotspots for Interaction with Alpha-1 Subunit:\n")
    for residue, role in interaction_hotspots.items():
        # Extract amino acid and position for clear output
        amino_acid = residue[0]
        position = residue[1:]
        print(f"Residue: {amino_acid}{position}")
        print(f"Role: {role}\n")

    print("-" * 70)

    # 2. Answer for gating modulation hotspots
    print("\n[Question 2] Hotspots for Fine-Tuning Gating Properties:\n")
    for residue, role in gating_modulation_hotspots.items():
        # Split residue string if it contains multiple residues
        residue_list = residue.split(', ')
        print_residue_str = ", ".join([f"{r[0]}{r[1:]}" for r in residue_list])
        print(f"Residue(s): {print_residue_str}")
        print(f"Role: {role}\n")


if __name__ == '__main__':
    find_calcium_channel_hotspots()
    # The final answer containing the specific residues is printed above.
    # For question 1: Y372, W391, E413, I415
    # For question 2: Y372, W391, E413, I415, W300, G302
    final_answer = {
        "interaction": ["Y372", "W391", "E413", "I415"],
        "gating_modulation": ["Y372", "W391", "E413", "I415", "W300", "G302"]
    }
    # This structure is just for internal representation and is fully displayed in the print output.
    # No direct return value is requested, only the printed output from the script.