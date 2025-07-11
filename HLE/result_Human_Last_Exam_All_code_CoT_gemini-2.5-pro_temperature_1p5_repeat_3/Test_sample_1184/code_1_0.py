import sys

def find_calcium_channel_hotspots():
    """
    This script identifies and prints key amino acid residues in the human
    voltage-gated calcium channel beta-1 subunit (CACNB1) for interaction
    and gating modulation of the alpha-1 subunit.
    """
    # Define the protein of interest for clarity
    protein_name = "Human Voltage-Gated Calcium Channel Beta-1 Subunit (CACNB1)"
    uniprot_id = "Q02641"

    # --- Data for Part 1: Interaction Hotspots ---
    # These residues form the core hydrophobic pocket in the Guanylate Kinase (GK)
    # domain that is essential for high-affinity binding to the alpha-1 subunit's
    # Alpha Interaction Domain (AID).
    interaction_hotspots = [
        {"residue": "Tyrosine", "abbreviation": "Tyr", "position": 421},
        {"residue": "Tryptophan", "abbreviation": "Trp", "position": 442},
        {"residue": "Isoleucine", "abbreviation": "Ile", "position": 444},
        {"residue": "Valine", "abbreviation": "Val", "position": 475},
        {"residue": "Tyrosine", "abbreviation": "Tyr", "position": 477},
    ]

    # --- Data for Part 2: Gating Modulation Hotspots ---
    # These residues are involved in the fine-tuning of channel kinetics
    # (activation and inactivation) and are located outside the main binding pocket.
    gating_hotspots = [
        {
            "residue": "KEKFEKQKE sequence",
            "abbreviation": "KEKFEKQKE",
            "position": "47-55",
            "note": "This N-terminal region is critical for modulating channel activation and inactivation.",
        },
        {
            "residue": "Aspartic acid",
            "abbreviation": "Asp",
            "position": 158,
            "note": "Located in the SH3 domain, this residue influences the voltage-dependence of channel activation.",
        },
    ]

    # --- Print the results ---
    print(f"Analysis for: {protein_name} (UniProt: {uniprot_id})\n")
    
    # Print results for Question 1
    print("=" * 75)
    print("1. Residues from Beta-1 Subunit for Interaction with Alpha-1 Subunit")
    print("=" * 75)
    print(
        "The following residues form the core binding pocket for the alpha-1 subunit:\n"
    )
    for res in interaction_hotspots:
        print(
            f"  - Residue: {res['residue']} ({res['abbreviation']}{res['position']})"
        )
    print("\n")

    # Print results for Question 2
    print("=" * 75)
    print("2. Residues from Beta-1 Subunit for Fine-Tuning Gating Properties")
    print("=" * 75)
    print(
        "The following residues and regions are critical for modulating channel gating:\n"
    )
    for res in gating_hotspots:
        print(
            f"  - Residue/Region: {res['residue']} ({res['abbreviation']}{res['position']})"
        )
        print(f"    - Note: {res['note']}")

if __name__ == '__main__':
    find_calcium_channel_hotspots()
    # Adding a marker for the final answer content as per instructions, although the
    # printed output above contains the complete answer. For a task like this,
    # the structured text output is the most valuable part of the answer.
    # The <<<>>> format is better suited for a single, concise value.
    # To comply, I will list the primary identified residues in the format.
    sys.stdout.write("<<<Interaction hotspots: Tyr421, Trp442, Ile444, Val475, Tyr477. Gating modulation hotspots: Region 47-55 (KEKFEKQKE), Asp158.>>>")