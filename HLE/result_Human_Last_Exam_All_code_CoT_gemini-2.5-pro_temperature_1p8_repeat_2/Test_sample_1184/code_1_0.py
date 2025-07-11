def print_channel_interaction_hotspots():
    """
    This script identifies and prints key amino acid residues on the human
    CaV beta-1 subunit (CACNB1) that are critical for its function.
    The information is based on structural and functional data from scientific literature.
    Residue numbering corresponds to human CACNB1, UniProt accession: P54283.
    """

    # --- Part 1: Residues for Interaction with alpha-1 Subunit ---
    # These residues are located in the Beta Interaction Domain (BID) and form the
    # binding groove for the alpha-1 subunit's Alpha Interaction Domain (AID).
    # This interaction is the primary anchor between the two subunits.
    interaction_hotspots = [
        ("Y", 216, "Forms hydrophobic contact in the binding pocket"),
        ("G", 218, "Backbone interaction"),
        ("W", 238, "Key hydrophobic contact"),
        ("V", 256, "Hydrophobic contact"),
        ("M", 261, "Hydrophobic contact"),
        ("R", 333, "Forms electrostatic interaction"),
        ("K", 336, "Forms electrostatic interaction"),
        ("D", 339, "Forms electrostatic interaction with AID"),
        ("L", 358, "Forms part of the hydrophobic groove"),
        ("F", 360, "Key hydrophobic contact"),
        ("R", 364, "Forms electrostatic interaction")
    ]

    print("---[ 1. Hotspots for Interaction with alpha-1 Subunit ]---\n")
    print("The following residues in the human beta-1 subunit's Beta Interaction Domain (BID)")
    print("are critical for high-affinity binding to the alpha-1 subunit:\n")
    print(f"{'Residue':<10} | {'Position':<10} | {'Note':<50}")
    print("-" * 75)
    for residue, pos, note in interaction_hotspots:
        print(f"{residue:<10} | {pos:<10} | {note:<50}")
    print("\n" + "="*75 + "\n")


    # --- Part 2: Residues for Fine-Tuning Gating Properties ---
    # Gating modulation is controlled by regions outside the core BID. These domains
    # influence the channel's opening, closing, and inactivation kinetics.
    gating_modulation_regions = {
        "N-terminal Domain": {
            "range": "1-20",
            "description": "This flexible N-terminus is a key modulator of voltage-dependent inactivation.",
            "example_residues": [("M", 1), ("S", 2), ("K", 3), ("S", 12), ("S", 20)]
        },
        "Hook Region": {
            "range": "approx. 166-189",
            "description": "This linker between the SH3 and GK domains strongly influences gating for CaV2.3 channels.",
            "example_residues": [("A", 166), ("E", 167), ("P", 168), ("P", 174), ("P", 180), ("D", 189)]
        }
    }

    print("---[ 2. Hotspots for Gating Modulation of alpha-1 Subunit ]---\n")
    print("Gating modulation is controlled by distinct functional regions rather than single scattered residues.")
    print("The primary hotspots for fine-tuning the alpha-1 subunit's gating properties are:\n")

    for region_name, details in gating_modulation_regions.items():
        print(f"Region: {region_name} (Residues {details['range']})")
        print(f"  Description: {details['description']}")
        print(f"  Example residues in this hotspot region include:")
        residue_str = ", ".join([f"{res}{pos}" for res, pos in details['example_residues']])
        print(f"    {residue_str}")
        print("-" * 75)


if __name__ == "__main__":
    print_channel_interaction_hotspots()
    print("\n<<<The output above contains the lists of requested residues and their positions.>>>")