def get_calcium_channel_hotspots():
    """
    This function identifies and prints key amino acid residues on the human
    voltage-gated calcium channel beta-1 subunit (CACNB1) involved in
    interaction and gating modulation with the alpha-1 subunit.

    The residue numbering is based on the canonical sequence from UniProt (ID: P54283).
    """

    # Part 1: Residues critical for binding to the alpha-1 subunit.
    # These residues are located in the Beta Interaction Domain (BID), a hydrophobic
    # groove that physically docks the Alpha Interaction Domain (AID) of the alpha-1 subunit.
    interaction_hotspots = {
        "Y244": "Forms a key part of the hydrophobic binding pocket for the AID.",
        "L248": "Contributes to the hydrophobic interactions essential for stable binding.",
        "V251": "A critical hydrophobic residue deep within the binding groove.",
        "L304": "Lines the wall of the binding groove, making contact with AID residues.",
        "V307": "Another key hydrophobic contact point for the AID.",
        "W368": "An essential, large hydrophobic residue that anchors the AID's Tryptophan (W).",
        "I372": "Contributes to the hydrophobic floor of the binding pocket."
    }

    # Part 2: Residues/regions important for fine-tuning the gating properties of the alpha-1 subunit.
    # These regions influence channel kinetics (e.g., activation, inactivation) allosterically.
    gating_modulation_hotspots = {
        "Region 31-71 (N-Terminus)": "This proximal N-terminal domain is a major determinant of the rate and voltage-dependence of channel inactivation.",
        "Region 170-180 (HOOK domain)": "Located in a variable loop, this region influences current density and inactivation kinetics, acting as a modulatory element.",
        "D253": "Located near the BID, this residue is implicated in transmitting the modulatory signal from the beta subunit to the alpha-1 voltage-sensing domains.",
        "E303": "Adjacent to a key binding residue (L304), this charged residue is thought to be involved in allosteric modulation of gating rather than just binding."
    }

    print("--- Hotspots on Human Beta-1 Subunit (CACNB1) for Alpha-1 Subunit Interaction ---")
    print("Based on UniProt Accession: P54283\n")

    print("1) Residues for Direct Interaction/Binding with Alpha-1 Subunit:")
    print("-" * 60)
    for residue, function in interaction_hotspots.items():
        # Extracting residue name and position for the final equation format
        res_name = residue[0]
        res_pos = residue[1:]
        print(f"Residue: {res_name} at Position: {res_pos} -> Role: {function}")
    print("-" * 60)

    print("\n2) Residues/Regions for Fine-Tuning Gating Properties of Alpha-1 Subunit:")
    print("-" * 60)
    for residue, function in gating_modulation_hotspots.items():
        if "Region" in residue:
            print(f"Region: {residue.split(' ')[1]} -> Role: {function}")
        else:
            res_name = residue[0]
            res_pos = residue[1:]
            print(f"Residue: {res_name} at Position: {res_pos} -> Role: {function}")
    print("-" * 60)
    print("\nNote: The final functional output of the channel is a result of complex allosteric interactions, and these represent key identified players.")

if __name__ == '__main__':
    get_calcium_channel_hotspots()