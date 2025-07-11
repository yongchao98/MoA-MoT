def get_calcium_channel_hotspots():
    """
    This script identifies and prints key amino acid residues in the human
    voltage-gated calcium channel beta-1 subunit (CACNB1, isoform beta-1b,
    UniProt P54283-2).

    It details residues critical for two distinct functions:
    1. Binding to the alpha-1 subunit.
    2. Fine-tuning the gating properties of the alpha-1 subunit.
    """

    # Data is based on UniProt entry P54283-2 (CACNB1, beta-1b isoform)
    # and a consensus from structural and functional studies.

    interaction_hotspots = [
        {'position': 313, 'aa_long': 'Tyrosine', 'aa_short': 'Y',
         'role': 'Forms a key part of the hydrophobic binding pocket for the alpha-1 subunit\'s AID motif.'},
        {'position': 324, 'aa_long': 'Tryptophan', 'aa_short': 'W',
         'role': 'Essential for high-affinity binding; its indole ring deeply inserts into a hydrophobic cleft on the alpha-1 subunit.'},
        {'position': 328, 'aa_long': 'Aspartate', 'aa_short': 'D',
         'role': 'Forms a critical hydrogen bond with a conserved Arginine in the alpha-1 subunit\'s AID, anchoring the interaction.'},
        {'position': 372, 'aa_long': 'Glutamate', 'aa_short': 'E',
         'role': 'Contributes to the electrostatic surface of the binding interface, stabilizing the complex.'},
        {'position': 376, 'aa_long': 'Valine', 'aa_short': 'V',
         'role': 'A conserved hydrophobic residue that contributes to the overall stability of the AID-BID complex.'}
    ]

    gating_modulation_hotspots = [
        {'position': '1-30', 'aa_long': 'N-terminal region', 'aa_short': 'N/A',
         'role': 'This entire unstructured region is a primary modulator of inactivation kinetics. Its presence accelerates the rate of channel inactivation.'},
        {'position': 16, 'aa_long': 'Arginine', 'aa_short': 'R',
         'role': 'A positively charged residue within the N-terminus. The charge distribution in this region is critical for its modulatory effect on gating.'},
        {'position': 18, 'aa_long': 'Arginine', 'aa_short': 'R',
         'role': 'Another key positively charged residue. Mutating these arginines can significantly alter voltage-dependent inactivation.'},
        {'position': 23, 'aa_long': 'Arginine', 'aa_short': 'R',
         'role': 'Part of a cluster of basic residues (R23, R24) that fine-tunes the interaction with other channel domains to control the gate.'}
    ]

    print("---[ Calcium Channel Beta-1 Subunit Analysis ]---\n")
    print("Residue numbering is based on human beta-1b isoform (UniProt: P54283-2).\n")

    print("1) Hotspots for Interaction with Alpha-1 Subunit:")
    print("-" * 50)
    print(f"{'Position':<10} {'Amino Acid':<15} {'Role'}")
    print(f"{'--------':<10} {'------------':<15} {'----'}")
    for residue in interaction_hotspots:
        aa_str = f"{residue['aa_long']} ({residue['aa_short']})"
        print(f"{residue['position']:<10} {aa_str:<15} {residue['role']}")
    print("-" * 50)

    print("\n2) Hotspots for Fine-Tuning Gating of Alpha-1 Subunit:")
    print("-" * 60)
    print(f"{'Position':<10} {'Amino Acid':<20} {'Role'}")
    print(f"{'--------':<10} {'------------------':<20} {'----'}")
    for residue in gating_modulation_hotspots:
        aa_str = f"{residue['aa_long']} ({residue['aa_short']})"
        print(f"{str(residue['position']):<10} {aa_str:<20} {residue['role']}")
    print("-" * 60)


if __name__ == '__main__':
    get_calcium_channel_hotspots()