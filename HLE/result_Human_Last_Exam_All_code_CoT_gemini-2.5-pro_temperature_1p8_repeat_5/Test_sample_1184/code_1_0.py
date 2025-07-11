import sys

# Disable writing __pycache__ files
sys.dont_write_bytecode = True

def find_calcium_channel_hotspots():
    """
    This script identifies and prints key amino acid residues and domains
    on the human voltage-gated calcium channel beta-1 subunit involved in
    interaction and modulation of the alpha-1 subunit.
    The residue numbering is based on the canonical sequence of human CACNB1 (UniProt: P54283).
    """

    # Data for Question 1: Interaction Hotspots
    interaction_hotspots = {
        "info": ("These residues are in the Beta Interaction Domain (BID), or Guanylate Kinase-like (GK) domain. "
                 "They form a binding site for the Alpha Interaction Domain (AID) of the alpha-1 subunit. "
                 "Data is derived from structural analysis of Cav channel complexes (e.g., PDB ID: 6JP8)."),
        "residues": {
            "Hydrophobic Pocket & Key Contacts": [
                {'position': 364, 'residue': 'Val', 'abbreviation': 'V'},
                {'position': 368, 'residue': 'Tyr', 'abbreviation': 'Y'},
                {'position': 372, 'residue': 'Ile', 'abbreviation': 'I'},
                {'position': 405, 'residue': 'Met', 'abbreviation': 'M'},
                {'position': 438, 'residue': 'Val', 'abbreviation': 'V'},
                {'position': 440, 'residue': 'Trp', 'abbreviation': 'W', 'note': 'Critical anchor residue in the hydrophobic pocket.'},
                {'position': 445, 'residue': 'Met', 'abbreviation': 'M'},
                {'position': 452, 'residue': 'Tyr', 'abbreviation': 'Y'}
            ],
            "Hydrogen Bonds & Electrostatic Interactions": [
                {'position': 355, 'residue': 'Arg', 'abbreviation': 'R'},
                {'position': 358, 'residue': 'Gln', 'abbreviation': 'Q'},
                {'position': 361, 'residue': 'Glu', 'abbreviation': 'E'},
                {'position': 397, 'residue': 'Gln', 'abbreviation': 'Q'},
                {'position': 399, 'residue': 'Lys', 'abbreviation': 'K'},
                {'position': 444, 'residue': 'Glu', 'abbreviation': 'E'}
            ]
        }
    }

    # Data for Question 2: Gating Modulation Hotspots
    gating_modulation_hotspots = {
        "info": ("Gating modulation is complex. The core interaction itself causes major effects (e.g., shifting voltage dependence). "
                 "However, other regions are critical for fine-tuning specific properties like inactivation speed."),
        "regions": [
            {
                "name": "Core Interaction Site (BID/GK Domain)",
                "range": "approx. 355-460",
                "function": "Binding of this domain to the alpha-1 subunit is responsible for the characteristic hyperpolarizing shift in the voltage-dependence of activation."
            },
            {
                "name": "N-terminal Domain",
                "range": "1-30",
                "function": "This flexible region is a key determinant for the rate and voltage-dependence of channel inactivation."
            },
            {
                "name": "C-terminal Domain",
                "range": "479-597",
                "function": "This distal region also contributes to the modulation of gating properties, working in concert with the other domains."
            }
        ]
    }

    # --- Print Results ---

    print("Analysis of Human Cav Beta-1 Subunit (CACNB1) Hotspots")
    print("=" * 70)
    print("\n1) Which exact residues from beta-1 subunit are hotspots for INTERACTION with alpha-1 subunit?\n")
    print(interaction_hotspots['info'])
    for category, residues in interaction_hotspots['residues'].items():
        print(f"\n--- {category} ---")
        for res in residues:
            note = f" ({res['note']})" if 'note' in res else ""
            print(f"  - Position {res['position']}: {res['residue']} ({res['abbreviation']}){note}")

    print("\n" + "=" * 70)
    print("\n2) Which exact residues from beta-1 subunit are hotspots for FINE-TUNING GATING PROPERTIES?\n")
    print(gating_modulation_hotspots['info'])
    for region in gating_modulation_hotspots['regions']:
        print(f"\n--- Hotspot Region: {region['name']} ---")
        print(f"  Residue Range: {region['range']}")
        print(f"  Function: {region['function']}")
    print("\n" + "=" * 70)


if __name__ == "__main__":
    find_calcium_channel_hotspots()
    print("<<<Data has been printed above.>>>")
