def get_calcium_channel_hotspots():
    """
    This script identifies and prints key amino acid residues on the human
    calcium channel beta-1 subunit (CACNB1) for its interaction with and
    modulation of the alpha-1 subunit (Cav2.3/CACNA1E).
    """

    # --- Data for Question 1 ---
    # These residues are located in the Beta Interaction Domain (BID), which is
    # part of the beta-1 subunit's core Guanylate Kinase (GK) domain.
    # They form a binding groove for the alpha-1 subunit's AID helix.
    # Residue numbering is based on the human CACNB1 canonical sequence (UniProt: P54289-1).
    interaction_hotspots = [
        ("Tyrosine", "Y", 211, "Forms a critical part of the hydrophobic binding pocket."),
        ("Lysine", "K", 247, "Forms a key electrostatic bond (salt bridge) with the AID."),
        ("Glutamic acid", "E", 251, "Contributes to electrostatic interactions with the AID helix."),
        ("Leucine", "L", 255, "Key hydrophobic interaction with a conserved residue in AID."),
        ("Isoleucine", "I", 259, "Contributes to the hydrophobic pocket that binds the AID."),
        ("Leucine", "L", 260, "Contributes to the hydrophobic pocket that binds the AID.")
    ]

    # --- Data for Question 2 ---
    # Gating modulation involves both the core interaction and other domains
    # that allosterically influence channel kinetics.
    gating_modulation_hotspots = [
        ("N-terminal Domain", "Region", "1-30", "Crucial for modulating the rate and voltage-dependence of inactivation."),
        ("GK Domain Residues", "e.g., Y211, K247, E251", "Binding at this core site itself modulates activation, typically shifting it to more negative voltages."),
        ("HOOK Region", "Region", "~155-184", "A linker region that contributes to isoform-specific differences in gating modulation.")
    ]

    # --- Print the results ---
    print("--- 1) Beta-1 Subunit Hotspots for Interaction with Alpha-1 Subunit ---")
    print("The key residues for direct binding are in the Beta Interaction Domain (BID):")
    print(f"{'Residue':<18} {'Position':<10} {'Function in Binding':<60}")
    print("-" * 90)
    for residue, code, pos, role in interaction_hotspots:
        res_str = f"{residue} ({code})"
        print(f"{res_str:<18} {pos:<10} {role:<60}")

    print("\n" * 2)

    print("--- 2) Beta-1 Subunit Hotspots for Fine-Tuning Gating Properties ---")
    print("Gating modulation is complex, involving several regions:")
    print(f"{'Domain/Region':<20} {'Position':<15} {'Role in Gating Modulation':<60}")
    print("-" * 95)
    for domain, pos, role_pos, role_gating in gating_modulation_hotspots:
        print(f"{domain:<20} {pos:<15} {role_gating:<60}")

if __name__ == '__main__':
    get_calcium_channel_hotspots()