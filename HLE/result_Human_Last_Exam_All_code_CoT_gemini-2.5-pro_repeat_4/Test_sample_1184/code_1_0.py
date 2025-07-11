import re

def print_residue_info(residue_dict, category_name):
    """
    Helper function to print residue information in a specific format.
    The position number will be printed digit by digit.
    """
    print(f"\n--- {category_name} ---\n")
    for residue, info in residue_dict.items():
        # Use regex to separate the amino acid letter from the position number
        match = re.match(r"([A-Z])(\d+)", residue)
        if not match:
            # Handle region-based info
            print(f"Region: {residue}\n  Role: {info['description']}\n")
            continue
            
        aa_letter = match.group(1)
        position_str = match.group(2)
        
        # Split the position number into individual digits
        position_digits = " ".join(list(position_str))
        
        aa_name = info['name']
        description = info['description']
        
        print(f"Residue: {aa_letter} ({aa_name}) at position {position_digits}")
        print(f"  Role: {description}\n")

def main():
    """
    Main function to define and print information about key residues in the
    human voltage-gated calcium channel beta-1 subunit (CACNB1).
    """
    # --- Data based on scientific literature and structural biology ---
    # UniProt ID for human CACNB1: P54289

    # 1) Residues critical for the physical interaction with the alpha-1 subunit.
    # These residues are in the Beta Interaction Domain (BID) within the GK domain.
    interaction_hotspots = {
        "Y216": {
            "name": "Tyrosine",
            "description": "Forms part of the conserved hydrophobic groove that docks the alpha-1 subunit's AID helix."
        },
        "E268": {
            "name": "Glutamic acid",
            "description": "Forms a critical salt bridge with a conserved Arginine residue in the alpha-1 AID, anchoring the interaction."
        },
        "W321": {
            "name": "Tryptophan",
            "description": "A key hydrophobic residue in the binding pocket, essential for high-affinity binding."
        },
        "I324": {
            "name": "Isoleucine",
            "description": "Contributes to the hydrophobic nature of the core binding site."
        },
        "F358": {
            "name": "Phenylalanine",
            "description": "Another crucial residue lining the hydrophobic binding groove for the alpha-1 AID."
        },
        "L361": {
            "name": "Leucine",
            "description": "Participates in the hydrophobic interactions that stabilize the alpha1-beta subunit complex."
        }
    }

    # 2) Residues/regions critical for fine-tuning the gating properties of the alpha-1 subunit.
    # These sites are often outside the primary binding interface and are involved in allosteric regulation.
    gating_modulation_hotspots = {
        "Residues 1-56": {
            "name": "N-Terminal Region",
            "description": "This variable region is a major determinant of the voltage-dependence and kinetics of channel inactivation."
        },
        "W184": {
            "name": "Tryptophan",
            "description": "Located in the SH3 domain. This residue is vital for the correct folding and structural integrity of the SH3 domain, which is required to properly translate beta subunit binding into modulation of channel activation."
        },
        "D236 & E237": {
            "name": "Aspartic acid & Glutamic acid",
            "description": "Part of an allosteric site on the GK domain (separate from the AID binding pocket) that influences isoform-specific gating modulation."
        },
        "R372 & K375": {
            "name": "Arginine & Lysine",
            "description": "Residues on the 'backside' of the GK domain, also implicated in transmitting the allosteric modulatory effects to the alpha-1 subunit."
        }
    }

    print("## Analysis of Human Beta-1 Subunit (CACNB1) Hotspot Residues ##\n")
    print("Based on structural data and mutagenesis studies, the following residues in the human")
    print("voltage-gated calcium channel beta-1 subunit (UniProt: P54289) are critical for its function.")

    print_residue_info(
        interaction_hotspots,
        "Question 1: Hotspots for Interaction with Alpha-1 Subunit"
    )
    print_residue_info(
        gating_modulation_hotspots,
        "Question 2: Hotspots for Fine-Tuning Gating Properties"
    )

if __name__ == "__main__":
    main()
<<<