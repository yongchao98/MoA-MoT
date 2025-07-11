def get_calcium_channel_hotspots():
    """
    This script identifies and prints key amino acid residues in the human
    voltage-gated calcium channel beta-1 subunit (CACNB1) responsible for
    interaction with and modulation of the alpha-1 subunit (Cav2.3).
    The residue numbering is based on the canonical human CACNB1 sequence (UniProt ID: Q02641).
    """

    # --- 1) Residues for Interaction with alpha-1 subunit ---
    # These residues form the hydrophobic binding groove in the Guanylate Kinase (GK) domain
    # that directly docks the Alpha Interaction Domain (AID) of the alpha-1 subunit.
    interaction_hotspots = [
        ("Leucine", "Leu", "L", 265),
        ("Tyrosine", "Tyr", "Y", 268),
        ("Valine", "Val", "V", 288),
        ("Isoleucine", "Ile", "I", 343),
        ("Tryptophan", "Trp", "W", 344),
        ("Methionine", "Met", "M", 368),
        ("Leucine", "Leu", "L", 371)
    ]

    # --- 2) Residues for Fine-Tuning Gating Properties ---
    # These residues and domains are critical for allosteric modulation of the
    # alpha-1 subunit's gating kinetics (e.g., activation, inactivation).
    # They are often located outside the primary binding interface.
    modulation_hotspots = [
        ("N-terminal Domain", "Region", "N/A", "1-45"),
        ("Tryptophan", "Trp", "W", 237)
    ]

    modulation_descriptions = {
        "1-45": "This entire domain is a key modulator of channel inactivation kinetics.",
        "237": "Located in the SH3 domain. Crucial for modulating voltage-dependence of inactivation, with less impact on direct binding affinity."
    }

    print("="*80)
    print("Hotspot Residues in Human Calcium Channel Beta-1 Subunit (CACNB1)")
    print("="*80)
    print("\n--- 1) Hotspots for Interaction with Alpha-1 Subunit (Binding) ---\n")
    print("These residues form the core binding pocket in the GK domain.")
    print(f"{'Amino Acid':<12} {'3-Letter':<10} {'1-Letter':<10} {'Position':<10}")
    print("-"*52)
    for name, code3, code1, pos in interaction_hotspots:
        print(f"{name:<12} {code3:<10} {code1:<10} {pos:<10}")

    print("\n\n--- 2) Hotspots for Fine-Tuning Gating Properties (Modulation) ---\n")
    print("These residues/domains are critical for allosteric regulation of channel function.")
    print(f"{'Residue/Domain':<20} {'Position':<10} {'Description'}")
    print("-"*80)
    for name, _, _, pos in modulation_hotspots:
        pos_str = str(pos)
        print(f"{name:<20} {pos_str:<10} {modulation_descriptions[pos_str]}")

if __name__ == '__main__':
    get_calcium_channel_hotspots()