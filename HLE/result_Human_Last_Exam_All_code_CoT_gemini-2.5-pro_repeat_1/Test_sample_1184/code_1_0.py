import textwrap

def get_calcium_channel_hotspots():
    """
    This script provides information on key residues in the human
    voltage-gated calcium channel beta-1 subunit (CACNB1) for its interaction
    with and modulation of the alpha-1 subunit (e.g., Cav2.3 / CACNA1E).

    The information is based on published structural and functional data.
    Residue numbering corresponds to the canonical isoform 1 (beta-1b) of
    human CACNB1 (UniProt ID: P54289-1).
    """

    # Question 1: Residues for interaction with the alpha-1 subunit
    # These residues are primarily in the Beta Interaction Domain (BID),
    # a conserved alpha-helix within the guanylate kinase-like (GK) domain
    # that forms the core binding surface with the alpha-1 subunit's AID.
    interaction_hotspots = {
        "K354 (Lysine 354)": "Forms a key electrostatic interaction (salt bridge) with the alpha-1 subunit's AID.",
        "Y357 (Tyrosine 357)": "Part of a conserved motif. Its aromatic ring makes critical hydrophobic contact with the AID groove.",
        "Y360 (Tyrosine 360)": "Like Y357, this tyrosine is crucial for the hydrophobic interaction anchoring the beta subunit to the alpha-1.",
        "L364 (Leucine 364)": "Contributes to the hydrophobic interface, stabilizing the core interaction.",
        "W368 (Tryptophan 368)": "Inserts deeply into a hydrophobic pocket on the AID, acting as a major anchor point.",
        "E371 (Glutamate 371)": "Forms an important salt bridge with a corresponding basic residue on the AID.",
        "V374 (Valine 374)": "Participates in the hydrophobic binding interface.",
        "I375 (Isoleucine 375)": "Another key hydrophobic residue that helps secure the beta subunit in the AID groove.",
        "K378 (Lysine 378)": "Forms a salt bridge at the periphery of the interaction site, contributing to binding affinity."
    }

    # Question 2: Residues/regions for fine-tuning gating properties
    # These regions modulate channel kinetics and voltage dependence, often
    # through allosteric mechanisms or interactions outside the primary BID-AID site.
    gating_modulation_hotspots = {
        "N-terminal Domain (residues 1-73)": "This entire domain is a major modulator. It is highly variable among splice forms and significantly influences the kinetics and voltage-dependence of channel inactivation. Its presence accelerates inactivation of Cav2.3 channels.",
        "P232 (Proline 232)": "Located in the flexible 'hinge' linker between the SH3 and GK domains. This residue is critical for the correct orientation between the two domains, which allosterically influences the channel's gating properties.",
        "W157 (Tryptophan 157)": "A conserved residue within the SH3 domain. While not at the binding interface, its mutation can disrupt the structural integrity of the SH3 domain, which in turn affects channel trafficking and can allosterically alter gating modulation."
    }

    print("="*80)
    print("Hotspot Residues in Human Calcium Channel Beta-1 Subunit (CACNB1)")
    print("="*80)
    print("\n1) Hotspots for Interaction with Alpha-1 Subunit (Binding):")
    print("-" * 60)
    for residue, description in interaction_hotspots.items():
        print(f"\nResidue: {residue}")
        print("  Role:")
        # Wrap the text for better readability
        wrapped_text = textwrap.fill(description, width=70, initial_indent='    ', subsequent_indent='    ')
        print(wrapped_text)

    print("\n\n2) Hotspots for Fine-Tuning Alpha-1 Subunit Gating (Modulation):")
    print("-" * 60)
    for region, description in gating_modulation_hotspots.items():
        print(f"\nRegion/Residue: {region}")
        print("  Role:")
        wrapped_text = textwrap.fill(description, width=70, initial_indent='    ', subsequent_indent='    ')
        print(wrapped_text)
    print("\n" + "="*80)

if __name__ == '__main__':
    get_calcium_channel_hotspots()