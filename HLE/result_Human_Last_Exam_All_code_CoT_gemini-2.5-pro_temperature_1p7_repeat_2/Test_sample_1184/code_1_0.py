def get_calcium_channel_hotspots():
    """
    This script identifies and prints key amino acid residues in the human
    voltage-gated calcium channel beta-1 subunit (CACNB1) for its interaction
    with and modulation of the alpha-1 subunit (CACNA1E).

    Residue numbering is based on the canonical human CACNB1 isoform 1 (UniProt: P54283).
    """

    # Part 1: Residues that are hotspots for the physical interaction with the alpha-1 subunit.
    # These residues are primarily in the Beta Interaction Domain (BID) and form a
    # hydrophobic groove that binds the Alpha Interaction Domain (AID) of the alpha-1 subunit.
    interaction_hotspots = {
        'W238': 'Tryptophan at position 238. Forms a key part of the hydrophobic binding pocket for the AID helix.',
        'Y240': 'Tyrosine at position 240. Contributes to the hydrophobic pocket and can form hydrogen bonds.',
        'E279': 'Glutamate at position 279. Forms electrostatic interactions (salt bridges) with the AID.',
        'I281': 'Isoleucine at position 281. A critical hydrophobic residue deep within the binding pocket.',
        'V285': 'Valine at position 285. Contributes to the hydrophobic nature of the binding pocket.',
        'L288': 'Leucine at position 288. Another key hydrophobic contact point for the AID helix.',
        'E330': 'Glutamate at position 330. Forms important salt bridges that stabilize the complex.'
    }

    # Part 2: Residues that are hotspots for fine-tuning the gating properties of the alpha-1 subunit.
    # This list includes the core interaction residues (as physical binding is required for modulation)
    # as well as residues in adjacent domains known to allosterically influence channel kinetics.
    gating_modulation_hotspots = {
        'Y214': 'Tyrosine at position 214. Located in the SH3 domain, mutations here affect channel inactivation kinetics.',
        'W238': 'Tryptophan at position 238. Disruption of binding at this site abolishes gating modulation.',
        'Y240': 'Tyrosine at position 240. Its interaction is crucial for transducing the modulatory signal.',
        'D251': 'Aspartate at position 251. Located in the "Hook" region connecting domains; affects voltage-dependent activation.',
        'I281': 'Isoleucine at position 281. Essential for the high-affinity interaction that underpins all gating modulation.',
        'E330': 'Glutamate at position 330. Proper electrostatic interaction here is required for normal modulation of gating.'
    }

    print("--- Part 1: Beta-1 Subunit Hotspots for INTERACTION with Alpha-1 Subunit ---")
    print("The following residues of the human beta-1 subunit (CACNB1) are critical for forming the high-affinity binding interface with the alpha-1 subunit:\n")
    for residue, description in interaction_hotspots.items():
        # Extracts the amino acid, number, and full name from the key and description
        aa_code = residue[0]
        position = residue[1:]
        print(f"Residue: {aa_code}{position} ({description.split(' at')[0]})")
        print(f"  - Position: {position}")
        print(f"  - Role: {description.split('. ')[1]}\n")


    print("\n--- Part 2: Beta-1 Subunit Hotspots for GATING MODULATION of Alpha-1 Subunit ---")
    print("The following residues of the human beta-1 subunit (CACNB1) are critical for fine-tuning the gating properties (e.g., activation, inactivation) of the alpha-1 subunit:\n")
    for residue, description in gating_modulation_hotspots.items():
        # Extracts the amino acid, number, and full name from the key and description
        aa_code = residue[0]
        position = residue[1:]
        print(f"Residue: {aa_code}{position} ({description.split(' at')[0]})")
        print(f"  - Position: {position}")
        print(f"  - Role: {description.split('. ')[1]}\n")

if __name__ == '__main__':
    get_calcium_channel_hotspots()