import textwrap

def print_channel_interaction_info():
    """
    Prints information about the key residues in the human voltage-gated calcium
    channel beta-1 subunit involved in interaction and gating modulation with the
    alpha-1 subunit.
    
    The information is based on structural biology and mutagenesis studies, with
    residue numbering corresponding to the human CACNB1 canonical sequence (UniProt P54283).
    """

    # 1) Residues for direct interaction with the alpha-1 subunit
    # These residues form the hydrophobic groove in the Beta Interaction Domain (BID)
    # on the GK domain of the beta-1 subunit, where the alpha-1 subunit's Alpha
    # Interaction Domain (AID) helix binds. The list is derived from the conserved
    # structure of the CaV alpha1/beta complex (e.g., PDB: 1T0J).
    binding_hotspots = {
        "Leu350 (L350)": "Forms part of the essential hydrophobic groove that docks the alpha-1 subunit.",
        "Val354 (V354)": "A key hydrophobic contact point within the binding groove.",
        "Tyr385 (Y385)": "Contributes to the binding pocket, interacting with the alpha-1 subunit's AID helix.",
        "Val393 (V393)": "Essential for the hydrophobic nature of the binding groove.",
        "Ile410 (I410)": "A critical residue deep within the groove, providing a snug fit for the alpha-1 AID. (Note: Often a Valine in other beta-subunits, but functionally equivalent).",
        "Leu415 (L415)": "Lines the wall of the binding groove, making direct hydrophobic contact.",
        "Leu418 (L418)": "Another key leucine residue that stabilizes the alpha-1/beta-1 interaction."
    }

    # 2) Residues for fine-tuning the gating properties of the alpha-1 subunit
    # These residues affect channel kinetics (activation/inactivation) without necessarily
    # disrupting the main binding interaction. They are identified from mutagenesis studies.
    gating_modulation_hotspots = {
        "Leu7 (L7)": "Located in the N-terminus. Mutation of this residue significantly alters the voltage-dependence of both activation and inactivation, demonstrating its role in fine-tuning channel gating (Takahashi et al., 2005, J. Physiol).",
        "Tyr100 (Y100)": "Located in the variable linker region (V1) between the N-terminus and the SH3 domain. Mutation of this residue was shown to dramatically alter inactivation kinetics, suggesting it is a key modulator of channel closing (Restituito et al., 2000, J. Biol. Chem)."
    }

    # -- Printing the Results --

    print("Analysis of Human CaV Beta-1 Subunit (CACNB1) Hotspot Residues")
    print("=" * 70)

    # Print Part 1: Interaction Hotspots
    print("\n1) Hotspot residues in the human beta-1 subunit for INTERACTION with the alpha-1 subunit:")
    print("-" * 70)
    for residue, description in binding_hotspots.items():
        print(f"Residue: {residue}")
        # textwrap is used for clean formatting of longer descriptions
        wrapped_desc = textwrap.fill(f"  Function: {description}", width=70, subsequent_indent='            ')
        print(wrapped_desc)
        print("-" * 20)

    # Print Part 2: Gating Modulation Hotspots
    print("\n2) Hotspot residues in the human beta-1 subunit for GATING MODULATION of the alpha-1 subunit:")
    print("-" * 70)
    for residue, description in gating_modulation_hotspots.items():
        print(f"Residue: {residue}")
        wrapped_desc = textwrap.fill(f"  Function: {description}", width=70, subsequent_indent='            ')
        print(wrapped_desc)
        print("-" * 20)

if __name__ == '__main__':
    print_channel_interaction_info()
