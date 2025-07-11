def get_calcium_channel_hotspots():
    """
    Identifies and prints key amino acid residues in the human voltage-gated
    calcium channel beta-1 subunit (CACNB1, UniProt: P54283).

    The script outlines two sets of hotspots:
    1. Residues crucial for the physical interaction with the alpha-1 subunit.
    2. Residues responsible for fine-tuning the gating properties of the channel.
    """

    # 1) Residues from beta-1 subunit that are hotspots for interaction with alpha-1 subunit
    # These residues are primarily located in the Beta Interaction Domain (BID),
    # a conserved region within the Guanylate Kinase (GK) domain of CACNB1.
    # They form a hydrophobic and electrostatic interface for the alpha-1 subunit's
    # Alpha Interaction Domain (AID).
    interaction_hotspots = [
        {'residue': 'W316', 'role': 'Forms deep hydrophobic pocket for AID binding'},
        {'residue': 'Y318', 'role': 'Contributes to the hydrophobic binding pocket'},
        {'residue': 'I357', 'role': 'Hydrophobic contact with the AID helix'},
        {'residue': 'F391', 'role': 'Key hydrophobic residue in the binding groove'},
        {'residue': 'Y393', 'role': 'Forms part of the hydrophobic wall of the binding pocket'},
        {'residue': 'W440', 'role': 'Critical hydrophobic anchor point for the AID helix'},
        {'residue': 'E441', 'role': 'Forms a key salt bridge with a conserved arginine on the AID'}
    ]

    # 2) Residues from beta-1 subunit that are hotspots for gating modulation
    # These residues or regions influence the channel's opening/closing kinetics
    # and voltage-dependence.
    modulation_hotspots = [
        {
            'residue': 'N-terminal Domain (residues 1-86)',
            'role': 'Significantly modulates the speed and voltage-dependence of channel inactivation'
        },
        {
            'residue': 'Y169',
            'role': 'Alters the voltage-dependence of channel activation. Mutation affects how readily the channel opens in response to voltage changes'
        },
        {
            'residue': 'D1-D2 Linker (residues 87-164)',
            'role': 'Contributes to the regulation of voltage-dependent inactivation kinetics'
        }
    ]

    print("--- 1. Hotspots for Interaction with alpha-1 Subunit ---")
    print("These residues in the beta-1 subunit's 'Beta Interaction Domain' (BID) are critical for binding to the alpha-1 subunit.\n")
    for spot in interaction_hotspots:
        print(f"Residue: {spot['residue']:<35} Role: {spot['role']}")

    print("\n\n--- 2. Hotspots for Gating Modulation of alpha-1 Subunit ---")
    print("These residues and domains in the beta-1 subunit fine-tune the gating behavior of the channel.\n")
    for spot in modulation_hotspots:
        print(f"Residue/Region: {spot['residue']:<35} Role: {spot['role']}")

# Execute the function to print the information
get_calcium_channel_hotspots()
<<<
--- 1. Hotspots for Interaction with alpha-1 Subunit ---
These residues in the beta-1 subunit's 'Beta Interaction Domain' (BID) are critical for binding to the alpha-1 subunit.

Residue: W316                                Role: Forms deep hydrophobic pocket for AID binding
Residue: Y318                                Role: Contributes to the hydrophobic binding pocket
Residue: I357                                Role: Hydrophobic contact with the AID helix
Residue: F391                                Role: Key hydrophobic residue in the binding groove
Residue: Y393                                Role: Forms part of the hydrophobic wall of the binding pocket
Residue: W440                                Role: Critical hydrophobic anchor point for the AID helix
Residue: E441                                Role: Forms a key salt bridge with a conserved arginine on the AID


--- 2. Hotspots for Gating Modulation of alpha-1 Subunit ---
These residues and domains in the beta-1 subunit fine-tune the gating behavior of the channel.

Residue/Region: N-terminal Domain (residues 1-86)   Role: Significantly modulates the speed and voltage-dependence of channel inactivation
Residue/Region: Y169                                Role: Alters the voltage-dependence of channel activation. Mutation affects how readily the channel opens in response to voltage changes
Residue/Region: D1-D2 Linker (residues 87-164)    Role: Contributes to the regulation of voltage-dependent inactivation kinetics
>>>