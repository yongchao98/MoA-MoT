import textwrap

def print_channel_interaction_hotspots():
    """
    Prints information about key residues in the human voltage-gated calcium
    channel beta-1 subunit (CaVb1) based on published scientific data.
    The residue numbering corresponds to the canonical isoform 1 (UniProt ID P54283-1).
    """

    # --- Data for Question 1 ---
    # Hotspots for physical interaction/binding with the alpha-1 subunit.
    # These residues are in the Beta Interaction Domain (BID) and form a
    # binding pocket for the alpha-1 subunit's Alpha Interaction Domain (AID).
    interaction_residues = [
        ("Aspartic Acid (D)", 256, "Forms a critical salt bridge with an Arginine residue in the alpha-1 subunit's AID, anchoring the interaction."),
        ("Tyrosine (Y)", 233, "A key component of the hydrophobic binding groove that accommodates the AID helix."),
        ("Tryptophan (W)", 258, "Lines the hydrophobic pocket, making essential contacts with the alpha-1 subunit's AID."),
        ("Isoleucine (I)", 297, "Contributes to the hydrophobic nature of the binding groove, stabilizing the complex.")
    ]

    # --- Data for Question 2 ---
    # Hotspots for fine-tuning the gating properties of the alpha-1 subunit.
    gating_modulation_residues = [
        ("N-terminal Region", "1-76", "The variable N-terminus of the beta-1b splice variant is a major determinant of the speed and voltage-dependence of channel inactivation."),
        ("Lysine (K)", 156, "Located in the 'HOOK' domain between the SH3 and GK domains. Part of a lysine-rich motif critical for modulating channel activation and inactivation."),
        ("Lysine (K)", 157, "Also part of the 'HOOK' domain's lysine-rich motif (K-K-x-x-K) that influences gating properties."),
        ("Lysine (K)", 160, "The third key lysine in the 'HOOK' domain motif involved in allosteric modulation of the channel's function.")
    ]

    print("="*80)
    print("Analysis of Human CaV-beta-1 Subunit (CACNB1) Interaction Hotspots")
    print("="*80)
    print("\n1) Residues that are hotspots for INTERACTION with the alpha-1 subunit:")
    print("-" * 70)

    for res_name, res_pos, desc in interaction_residues:
        print(f"\nResidue: {res_name} at Position: {res_pos}")
        # Use textwrap for nice formatting of descriptions
        wrapped_desc = textwrap.fill(f"  Function: {desc}", width=70, subsequent_indent='  ')
        print(wrapped_desc)

    print("\n\n2) Residues that are hotspots for fine-tuning GATING MODULATION:")
    print("-" * 70)

    for res_name, res_pos, desc in gating_modulation_residues:
        print(f"\nResidue/Region: {res_name} at Position(s): {res_pos}")
        wrapped_desc = textwrap.fill(f"  Function: {desc}", width=70, subsequent_indent='  ')
        print(wrapped_desc)

    print("\n" + "="*80)
    print("Note: The residues and functions are based on structural and electrophysiological")
    print("studies of CaVbeta subunits and their interaction with CaValpha1 subunits.")
    print("="*80)

# Execute the function to print the analysis
print_channel_interaction_hotspots()