import sys

def find_calcium_channel_hotspots():
    """
    This script identifies and prints key amino acid residues in the human
    CaV beta-1 subunit (CACNB1) based on established scientific literature
    and structural biology data.

    The residue numbering corresponds to the canonical human sequence (UniProt: P54283).
    """

    print("Identifying key residues in Human Voltage-Gated Calcium Channel Beta-1 Subunit (CACNB1, UniProt: P54283).")
    print("-" * 80)

    # 1. Residues that are hotspots for interaction with the alpha-1 subunit
    # These residues form the core of the AID-binding pocket in the GK domain.
    interaction_hotspots = [
        ('Y', 247, "Forms a key part of the hydrophobic binding pocket for the alpha-1 AID."),
        ('L', 251, "Contributes to the hydrophobic character of the AID binding pocket."),
        ('V', 312, "A critical residue lining the AID binding pocket."),
        ('L', 314, "Another key hydrophobic residue for high-affinity AID binding."),
        ('W', 372, "A large, essential residue that deeply anchors the alpha-1 AID helix."),
        ('M', 401, "Contributes to the hydrophobic surface of the binding groove."),
        ('I', 404, "Lines the binding pocket and makes critical hydrophobic contacts with the AID.")
    ]

    print("\n1) Hotspots for interaction with the alpha-1 subunit (AID-Binding Pocket):\n")
    for residue, position, description in interaction_hotspots:
        print(f"Residue: {residue}{position}\n  - Role: {description}\n")

    # 2. Residues that are hotspots for fine-tuning gating properties
    # These residues are involved in allosteric regulation of the alpha-1 subunit.
    gating_hotspots = [
        ('K', 20, "Part of a dibasic K20-R21 motif in the N-terminus, critical for slowing channel inactivation."),
        ('R', 21, "Part of a dibasic K20-R21 motif in the N-terminus, critical for slowing channel inactivation."),
        ('E', 192, "Located in the SH3 domain, forms a key salt bridge with R430 in the GK domain."),
        ('R', 430, "Located in the GK domain, its interaction with E192 maintains the beta subunit's conformation, which is essential for proper gating modulation.")
    ]

    print("-" * 80)
    print("\n2) Hotspots for fine-tuning gating properties of the alpha-1 subunit:\n")
    for residue, position, description in gating_hotspots:
        print(f"Residue: {residue}{position}\n  - Role: {description}\n")

if __name__ == '__main__':
    find_calcium_channel_hotspots()