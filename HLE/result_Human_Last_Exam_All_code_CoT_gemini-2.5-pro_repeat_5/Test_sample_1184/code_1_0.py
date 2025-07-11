import collections

def display_channel_residues():
    """
    Identifies and displays key residues on the human CaV beta-1 subunit (CACNB1)
    for interaction with and modulation of the CaV alpha-1 subunit.
    Residue numbering is based on the canonical human CACNB1 sequence (UniProt P54289, 597 aa).
    """

    # 1) Residues for direct interaction with the alpha-1 subunit
    # These residues form the hydrophobic binding groove in the Guanylate Kinase (GK) domain
    # that accommodates the Alpha Interaction Domain (AID) from the alpha-1 subunit.
    interaction_hotspots = collections.OrderedDict([
        (224, ("Leu", "L")),
        (226, ("Tyr", "Y")),
        (273, ("Val", "V")),
        (275, ("Ala", "A")),
        (277, ("Ile", "I")),
        (313, ("Leu", "L")),
        (316, ("Trp", "W")),
        (320, ("Trp", "W")),
        (340, ("Met", "M")),
        (344, ("Leu", "L")),
        (373, ("Tyr", "Y")),
        (374, ("Ala", "A")),
        (377, ("Ile", "I")),
    ])

    # 2) Residues for fine-tuning the gating properties of the alpha-1 subunit
    # This is a more complex process involving several regions.

    # A) A key motif in the N-terminus of the beta-1b splice variant acts as a "gating brake",
    # slowing channel activation.
    gating_brake_motif = collections.OrderedDict([
        (3, ("Lys", "K")),
        (4, ("Arg", "R")),
        (5, ("Arg", "R")),
        (7, ("Lys", "K")),
        (8, ("Lys", "K")),
        (9, ("Arg", "R")),
    ])

    # B) The "Hook" region of the GK domain interacts with the I-II loop of the alpha-1 subunit
    # (outside the main AID helix) and is critical for modulating voltage-dependent inactivation.
    hook_region = collections.OrderedDict([
        (351, ("Glu", "E")),
        (352, ("Asp", "D")),
        (355, ("Arg", "R")),
    ])

    print("---")
    print("1) Interaction Hotspots on Human Beta-1 Subunit (CACNB1)")
    print("These residues form the hydrophobic groove in the Guanylate Kinase (GK) domain that binds the Alpha Interaction Domain (AID) of the alpha-1 subunit.")
    print("-" * 20)
    for pos, (name, code) in interaction_hotspots.items():
        print(f"Residue: {name} ({code}) at position {pos}")

    print("\n" + "---")
    print("2) Gating Modulation Hotspots on Human Beta-1 Subunit (CACNB1)")
    print("Modulation is complex, but specific residues have been identified in key regions:")

    print("\nA) N-terminal 'Gating Brake' Motif (identified in splice variant Beta-1b):")
    print("These residues form a polybasic motif that slows down activation kinetics.")
    print("-" * 20)
    for pos, (name, code) in gating_brake_motif.items():
        print(f"Residue: {name} ({code}) at position {pos}")

    print("\nB) 'Hook' Region of the GK Domain:")
    print("This region interacts with the alpha-1 I-II loop and influences voltage-dependent inactivation.")
    print("-" * 20)
    for pos, (name, code) in hook_region.items():
        print(f"Residue: {name} ({code}) at position {pos}")
    print("---")


if __name__ == "__main__":
    display_channel_residues()
