import textwrap

def print_wrapped(text):
    """Helper function to print text with wrapping."""
    print(textwrap.fill(text, width=80))
    print("-" * 80)

def find_cav_beta_subunit_hotspots():
    """
    Identifies and prints key residues in the human calcium channel beta-1 subunit
    (CACNB1) for interaction with and modulation of the alpha-1 subunit.

    The information is compiled from conserved structural and functional data
    from sources like UniProt (P54283 for human CACNB1), PDB, and scientific literature.
    """
    print("*" * 80)
    print("Identifying Residue Hotspots in Human Voltage-Gated Calcium Channel Beta-1 Subunit (CACNB1)")
    print("*" * 80)

    # --- 1) Hotspots for Interaction with Alpha-1 Subunit ---
    print("\n[1] Hotspots for Interaction with the Alpha-1 Subunit")
    print("-" * 80)
    print_wrapped(
        "These residues are located in the Guanylate Kinase (GK) domain of the "
        "beta-1 subunit. They form a deep, conserved hydrophobic groove "
        "(often called the Beta Interaction Domain or BID) that physically "
        "binds the Alpha Interaction Domain (AID) of the alpha-1 subunit. "
        "This binding is essential for the formation and trafficking of the channel complex."
    )

    interaction_residues = {
        "Tyrosine-213 (Tyr, Y213)": "Forms a key part of the binding pocket wall.",
        "Glutamate-215 (Glu, E215)": "Participates in hydrogen bonding at the groove's rim.",
        "Leucine-217 (Leu, L217)": "Contributes to the hydrophobic interactions within the pocket.",
        "Tryptophan-238 (Trp, W238)": "A critical residue at the base of the binding pocket, making extensive contact with the AID helix.",
        "Valine-259 (Val, V259)": "Lines the hydrophobic groove, interacting with the alpha-helix of the AID.",
        "Phenylalanine-260 (Phe, F260)": "Another key hydrophobic residue lining the binding pocket.",
        "Leucine-297 (Leu, L297)": "Contributes to the hydrophobic surface of the pocket.",
        "Isoleucine-301 (Ile, I301)": "Lines the pocket, providing a surface for hydrophobic interaction."
    }

    print("The exact hotspot residues from the beta-1 subunit are:")
    for residue, description in interaction_residues.items():
        print(f"  - {residue}: {description}")

    # --- 2) Hotspots for Gating Modulation of Alpha-1 Subunit ---
    print("\n[2] Hotspots for Fine-Tuning Gating Properties of the Alpha-1 Subunit")
    print("-" * 80)
    print_wrapped(
        "Gating modulation is a complex process. It is fundamentally dependent on the "
        "stable interaction described above, so the binding hotspot residues are also "
        "critical for modulation. However, specific residues outside the immediate "
        "binding pocket have been shown to allosterically fine-tune gating kinetics, "
        "such as the speed and voltage-dependence of activation and inactivation."
    )

    modulation_residues = {
        "All residues listed in Part [1]": "Proper binding is the primary requirement for all gating modulation.",
        "Aspartate-152 (Asp, D152)": "Located in the SH3 domain. This residue and its neighbors are involved in an intramolecular interaction with the GK domain that is critical for modulating voltage-dependent inactivation.",
        "Lysine-154 (Lys, K154)": "Also in the SH3 domain, this residue is part of a key surface that influences communication between the SH3 and GK domains, thereby affecting channel inactivation.",
        "The SH3-GK Linker (residues ~165-205)": "While not a single residue, the flexibility and conformation of this entire linker region are crucial for transmitting the modulatory effects of the beta subunit to the alpha-1 subunit."
    }

    print("The exact hotspot residues and regions for gating modulation are:")
    for residue, description in modulation_residues.items():
        print(f"  - {residue}: {description}")
    
    print("\n" + "*" * 80)


if __name__ == '__main__':
    find_cav_beta_subunit_hotspots()