def find_calcium_channel_hotspots():
    """
    This script identifies and prints key amino acid residues on the human
    voltage-gated calcium channel beta-1 subunit (CACNB1, UniProt: P54283).
    The residues are categorized into two groups:
    1. Hotspots for the physical interaction with the alpha-1 subunit.
    2. Hotspots for the fine-tuning of the alpha-1 subunit's gating properties.
    The residue numbering is based on the canonical sequence of human CACNB1.
    """

    # 1) Residues critical for the physical interaction with the alpha-1 subunit.
    # These residues are located in the Beta Interaction Domain (BID) and form a
    # hydrophobic pocket that binds the Alpha Interaction Domain (AID) of the alpha-1 subunit.
    # Mutating these residues typically disrupts or abolishes the binding.
    interaction_hotspots = {
        "Y259": "Tyrosine at position 259",
        "W313": "Tryptophan at position 313",
        "L315": "Leucine at position 315",
        "I317": "Isoleucine at position 317",
        "V349": "Valine at position 349",
        "Y352": "Tyrosine at position 352"
    }

    # 2) Residues/domains critical for modulating the gating properties of the alpha-1 subunit.
    # These regions influence the channel's opening/closing kinetics. This is an
    # allosteric effect, distinct from the primary binding itself.
    gating_modulation_hotspots = {
        "Region 1-48": "The N-terminal domain (specifically in isoform CACNB1b) is a key determinant of channel inactivation speed.",
        "R108": "Arginine at position 108 (in the SH3 domain)",
        "K109": "Lysine at position 109 (in the SH3 domain)",
        "E110": "Glutamate at position 110 (in the SH3 domain)"
    }

    print("--- Human Calcium Channel Beta-1 Subunit (CACNB1) Analysis ---")

    print("\n[1] Hotspot Residues for Interaction with Alpha-1 Subunit:")
    print("----------------------------------------------------------")
    print("These residues form the core of the binding pocket for the alpha-1 subunit.")
    for residue, description in interaction_hotspots.items():
        print(f"- {residue}: {description}")

    print("\n[2] Hotspot Residues for Gating Modulation of Alpha-1 Subunit:")
    print("----------------------------------------------------------------")
    print("These residues and regions are crucial for fine-tuning the channel's function.")
    for residue, description in gating_modulation_hotspots.items():
        print(f"- {residue}: {description}")


if __name__ == "__main__":
    find_calcium_channel_hotspots()