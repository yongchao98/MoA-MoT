def find_hotspot_residues():
    """
    This function identifies and prints key amino acid residues on the human
    voltage-gated calcium channel beta-1 subunit (CACNB1, UniProt: P54289)
    for interaction with the alpha-1 subunit and for modulating its gating properties.
    """

    # --- 1) Interaction Hotspots with alpha-1 subunit ---
    # These residues form the core of the Beta Interaction Domain (BID), a hydrophobic
    # groove within the Guanylate Kinase-like (GK) domain. This groove directly
    # binds the Alpha Interaction Domain (AID) from the alpha-1 subunit.
    interaction_hotspots = {
        "Y249": "Tyrosine at position 249",
        "L252": "Leucine at position 252",
        "E256": "Glutamate at position 256",
        "V259": "Valine at position 259",
        "W335": "Tryptophan at position 335",
        "D336": "Aspartate at position 336",
        "E338": "Glutamate at position 338",
        "I340": "Isoleucine at position 340",
    }

    # --- 2) Gating Modulation Hotspots for alpha-1 subunit ---
    # These residues are not in the primary binding interface but have been shown
    # to allosterically influence the channel's gating kinetics (e.g., voltage
    # dependence of activation or inactivation).
    gating_modulation_hotspots = {
        "SH3-GK Linker (approx. 231-234)": "This flexible region, including Proline 230 (P230), connects two key domains and its conformation is critical for gating modulation.",
        "R350": "Arginine at position 350 - Influences the voltage-dependence of channel activation.",
        "E409": "Glutamate at position 409 - Affects the voltage-dependence of channel inactivation.",
        "R412": "Arginine at position 412 - Also affects the voltage-dependence of channel inactivation.",
    }

    print("--- Analysis of Human Calcium Channel Beta-1 Subunit (CACNB1) ---")
    print("\n")

    print("1. Hotspot Residues for Interaction with Alpha-1 Subunit:")
    print("----------------------------------------------------------")
    for residue, description in interaction_hotspots.items():
        aa = residue[0]
        pos = residue[1:]
        print(f"Residue: {aa} (Amino Acid) at Position: {pos} - ({description})")

    print("\n")
    print("2. Hotspot Residues for Gating Modulation of Alpha-1 Subunit:")
    print("---------------------------------------------------------------")
    for residue, description in gating_modulation_hotspots.items():
        if " " in residue: # Handle the conceptual region
            print(f"Region: {residue} - {description}")
        else:
            aa = residue[0]
            pos = residue[1:]
            print(f"Residue: {aa} (Amino Acid) at Position: {pos} - ({description})")

if __name__ == '__main__':
    find_hotspot_residues()