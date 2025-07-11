def get_calcium_channel_hotspots():
    """
    This function provides information on key amino acid residues in the human
    voltage-gated calcium channel beta-1 subunit (CACNB1).

    The residue numbers are based on the canonical sequence from UniProt (Accession: Q02641).
    """

    # --- Question 1: Interaction Hotspots with Alpha-1 Subunit ---
    # These residues are primarily located in the Guanylate Kinase (GK) domain and
    # form the hydrophobic binding pocket for the Alpha Interaction Domain (AID)
    # of the alpha-1 subunit. Data is derived from structural studies.
    interaction_hotspots = {
        "Description": "Key residues in the beta-1 subunit (CACNB1) forming the binding pocket for the alpha-1 subunit.",
        "Residues": {
            "L209": "Lines the binding pocket.",
            "Y213": "Forms a key hydrogen bond and hydrophobic interactions.",
            "V260": "Contributes to the hydrophobic pocket.",
            "I261": "Contributes to the hydrophobic pocket.",
            "W266": "Makes critical hydrophobic contact with the AID helix.",
            "M286": "Central residue in the hydrophobic pocket.",
            "L289": "Lines the binding pocket.",
            "I292": "Makes critical hydrophobic contact with the AID helix.",
            "Y325": "Contributes to the binding interface.",
            "W345": "Makes critical hydrophobic contact with the AID helix.",
            "L348": "Lines the binding pocket.",
            "F350": "Contributes to the hydrophobic pocket.",
            "I419": "Lines the binding pocket.",
            "I423": "Contributes to the binding interface."
        }
    }

    # --- Question 2: Gating Modulation Hotspots ---
    # These residues are involved in fine-tuning the gating properties (e.g., activation,
    # inactivation) of the alpha-1 subunit. They are often located outside the
    # primary binding pocket and act via allosteric mechanisms.
    gating_modulation_hotspots = {
        "Description": "Key residues in the beta-1 subunit (CACNB1) for fine-tuning the gating properties of the alpha-1 subunit.",
        "Residues": {
            # N-terminus residues
            "L9": "Located in the N-terminus; important for modulating current density and inactivation kinetics.",
            "M10": "Located in the N-terminus; important for modulating current density and inactivation kinetics.",
            # SH3-GK interface residues for allosteric control
            "K171": "In the SH3 domain; part of the SH3-GK interface crucial for voltage-dependent activation.",
            "D172": "In the SH3 domain; part of the SH3-GK interface crucial for voltage-dependent activation.",
            "K482": "In the GK domain; part of the SH3-GK interface crucial for voltage-dependent activation."
        }
    }

    # --- Print the results ---
    print("Residue analysis for Human Voltage-Gated Calcium Channel Beta-1 Subunit (CACNB1, UniProt: Q02641)\n")

    print("="*80)
    print("1) Hotspots for Interaction with Alpha-1 Subunit")
    print("="*80)
    print(interaction_hotspots["Description"])
    print("-" * 50)
    for residue, function in interaction_hotspots["Residues"].items():
        # Using string formatting to output the "equation-like" format
        print(f"Residue: {residue:<5s} | Function: {function}")
    print("\n")


    print("="*80)
    print("2) Hotspots for Fine-Tuning Gating Properties of Alpha-1 Subunit")
    print("="*80)
    print(gating_modulation_hotspots["Description"])
    print("-" * 50)
    for residue, function in gating_modulation_hotspots["Residues"].items():
        print(f"Residue: {residue:<5s} | Function: {function}")

if __name__ == '__main__':
    get_calcium_channel_hotspots()
    # Combine the answer into a single string for the final output as requested
    final_answer_string = (
        "Interaction Hotspots (Beta-1 Subunit): L209, Y213, V260, I261, W266, M286, L289, I292, Y325, W345, L348, F350, I419, I423. "
        "Gating Modulation Hotspots (Beta-1 Subunit): L9, M10, K171, D172, K482."
    )
    # The final answer format as requested.
    print(f"\n<<<Interaction Hotspots (Beta-1 Subunit): L209, Y213, V260, I261, W266, M286, L289, I292, Y325, W345, L348, F350, I419, I423. Gating Modulation Hotspots (Beta-1 Subunit): L9, M10, K171, D172, K482.>>>")
