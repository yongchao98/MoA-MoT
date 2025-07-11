def find_channel_interaction_residues():
    """
    This script identifies and prints key amino acid residues on the human
    CaV beta-1 subunit for interaction and gating modulation of the
    CaV2.3 alpha-1 subunit.

    Residue numbering is based on the canonical human CACNB1 isoform (UniProt: P54283).
    """

    # --- 1. Residues for physical interaction with the alpha-1 subunit ---
    # These residues are in the Beta Interaction Domain (BID), a hydrophobic pocket
    # on the GK domain that binds the Alpha Interaction Domain (AID) of the alpha-1 subunit.
    interaction_hotspots = {
        'Description': "Hotspots for direct binding to the alpha-1 subunit (AID-BID interaction)",
        'Residues': [
            {'residue': 'Y237', 'role': 'Forms wall of the binding pocket'},
            {'residue': 'W284', 'role': 'Key hydrophobic anchor for the AID helix'},
            {'residue': 'M313', 'role': 'Contributes to the hydrophobic pocket'},
            {'residue': 'F316', 'role': 'Contributes to the hydrophobic pocket'},
            {'residue': 'I317', 'role': 'Contributes to the hydrophobic pocket'},
            {'residue': 'Y370', 'role': 'Forms wall of the binding pocket'},
            {'residue': 'W373', 'role': 'Key hydrophobic anchor for the AID helix'},
            {'residue': 'M395', 'role': 'Contributes to the hydrophobic pocket'},
            {'residue': 'A399', 'role': 'Contributes to the hydrophobic pocket'},
        ]
    }

    # --- 2. Residues for fine-tuning gating properties of the alpha-1 subunit ---
    # Gating modulation is complex. Besides the BID, other regions allosterically
    # influence channel function. Here are two well-documented sites.
    gating_modulation_hotspots = [
        {
            'Region': "SH3 Domain Surface",
            'Description': "This cluster of basic residues on the SH3 domain surface is involved in modulating the kinetics of voltage-dependent inactivation.",
            'Residues': [
                {'residue': 'K181', 'role': 'Surface charge modulation'},
                {'residue': 'R183', 'role': 'Surface charge modulation'},
                {'residue': 'K185', 'role': 'Surface charge modulation'},
            ]
        },
        {
            'Region': "GK Domain Intramolecular 'beta-ANCHOR' Site",
            'Description': "These residues anchor the N-terminus of the beta-1 subunit to its own GK domain. This intramolecular interaction is critical for the specific modulatory effects of the beta-1 subunit on channel gating.",
            'Residues': [
                {'residue': 'K377', 'role': 'N-terminus anchor point'},
                {'residue': 'R469', 'role': 'N-terminus anchor point'},
            ]
        }
    ]

    # --- Print the results ---
    print("<<<")
    print("Based on published structural and functional data for the human CaV beta-1 subunit (UniProt: P54283):\n")

    # Print Interaction Hotspots
    print("1) Residues from beta-1 subunit for INTERACTION with alpha-1 subunit:")
    print("-" * 70)
    print(f"Location: Beta Interaction Domain (BID) on the GK-like domain")
    print(f"Description: {interaction_hotspots['Description']}.")
    print("Residues:")
    for r in interaction_hotspots['Residues']:
        print(f"  - {r['residue']}: {r['role']}")
    print("\n")

    # Print Gating Modulation Hotspots
    print("2) Residues from beta-1 subunit for GATING MODULATION of alpha-1 subunit:")
    print("-" * 70)
    print("Gating modulation is controlled by multiple regions. Key sites for fine-tuning include:\n")
    for area in gating_modulation_hotspots:
        print(f"Location: {area['Region']}")
        print(f"Description: {area['Description']}")
        print("Residues:")
        for r in area['Residues']:
            print(f"  - {r['residue']}: {r['role']}")
        print("") # Add a newline for spacing

    # Conclude with a note that interaction residues also affect modulation
    print("Note: The interaction hotspots listed in (1) are also fundamentally required for gating modulation, as binding is a prerequisite for any modulatory effect.")
    print(">>>")

if __name__ == '__main__':
    find_channel_interaction_residues()