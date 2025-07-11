def get_calcium_channel_hotspots():
    """
    This script provides information on the key residues of the human voltage-gated
    calcium channel beta-1 subunit involved in interacting with and modulating the
    alpha-1 subunit. The information is based on published structural and
    functional studies. Residue numbering is based on the human beta-1b isoform
    (UniProt: P54283-2).
    """

    # --- 1) Hotspots for Interaction (Binding) with alpha-1 subunit ---
    # These residues form a highly conserved hydrophobic pocket in the Guanylate
    # Kinase (GK) domain of the beta-1 subunit. This pocket is the primary
    # binding site for the Alpha Interaction Domain (AID) of the alpha-1 subunit.
    binding_hotspots = [
        {'residue': 'Y (Tyrosine)', 'position': 245, 'role': 'Forms a key part of the hydrophobic binding groove for the AID.'},
        {'residue': 'L (Leucine)', 'position': 251, 'role': 'Contributes to the hydrophobic surface of the binding pocket.'},
        {'residue': 'L (Leucine)', 'position': 340, 'role': 'A critical hydrophobic contact point for the AID helix.'},
        {'residue': 'I (Isoleucine)', 'position': 342, 'role': 'Contributes to the hydrophobic surface of the binding pocket.'},
        {'residue': 'W (Tryptophan)', 'position': 362, 'role': 'Provides a large hydrophobic surface for a deep and stable interaction.'},
        {'residue': 'L (Leucine)', 'position': 369, 'role': 'A critical hydrophobic contact point for the AID helix.'}
    ]

    print("1) Hotspots on Beta-1 Subunit for Interaction (Binding) with Alpha-1 Subunit:\n")
    for spot in binding_hotspots:
        print(f"Residue: {spot['residue']} at Position: {spot['position']}")
        print(f"Role: {spot['role']}\n")
    print("-" * 50)

    # --- 2) Hotspots for Fine-Tuning Gating Properties of alpha-1 subunit ---
    # Gating modulation is more complex and is often attributed to entire domains
    # rather than single residues. These domains allosterically influence the
    # alpha-1 subunit's voltage sensors and pore after the beta subunit has bound.
    gating_hotspots = [
        {'region': 'N-terminal Domain', 'positions': 'approx. 1-30',
         'role': 'This variable domain is crucial for modulating the kinetics and voltage-dependence of channel inactivation.'},
        {'region': 'HOOK Domain', 'positions': 'approx. 260-315',
         'role': 'This linker region between the SH3 and GK domains is a key determinant for modulating the voltage-dependence of channel activation.'}
    ]

    print("\n2) Hotspots on Beta-1 Subunit for Gating Modulation of Alpha-1 Subunit:\n")
    for spot in gating_hotspots:
        print(f"Region: {spot['region']} at Positions: {spot['positions']}")
        print(f"Role: {spot['role']}\n")

    print("\nNote: The information is compiled from scientific literature. 'Binding' refers to the physical attachment, while 'Gating Modulation' refers to the functional effect on channel opening and closing after attachment.")

if __name__ == '__main__':
    get_calcium_channel_hotspots()