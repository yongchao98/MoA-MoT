def get_calcium_channel_hotspots():
    """
    This function defines and prints the critical residues in the human
    calcium channel beta-1 subunit for interaction and gating modulation
    of the alpha-1 subunit.

    The information is compiled from established biophysical and structural biology studies.
    Residue numbering corresponds to UniProt accession P54289 (CACNB1_HUMAN).
    """

    interaction_hotspots = {
        'Y253 (Tyrosine)': 'Located in the GK domain; forms key interactions in the hydrophobic binding groove for the alpha-1 AID.',
        'M275 (Methionine)': 'Contributes to the hydrophobic binding pocket in the GK domain.',
        'V278 (Valine)': 'Contributes to the hydrophobic binding pocket in the GK domain.',
        'W282 (Tryptophan)': 'A critical residue in the GK domain, forms extensive hydrophobic contacts with the alpha-1 AID helix.',
        'L312 (Leucine)': 'Contributes to the hydrophobic binding pocket in the GK domain.',
        'A315 (Alanine)': 'Contributes to the hydrophobic binding pocket in the GK domain.',
        'L319 (Leucine)': 'Contributes to the hydrophobic binding pocket in the GK domain.',
        'I320 (Isoleucine)': 'Contributes to the hydrophobic binding pocket in the GK domain.',
        'V323 (Valine)': 'Contributes to the hydrophobic binding pocket in the GK domain.'
    }

    gating_modulation_hotspots = {
        'E134 (Glutamic acid)': 'Located on the surface of the SH3 domain; part of a proposed secondary site influencing gating modulation.',
        'E136 (Glutamic acid)': 'Located on the surface of the SH3 domain, adjacent to E134, contributing to gating modulation.',
        'W161 (Tryptophan)': 'Located in the SH3 domain. Mutation of this residue can uncouple gating modulation from high-affinity binding, highlighting its specific role in regulation.',
        'HOOK Region (~219-245)': 'This flexible linker between the SH3 and GK domains is not a single hotspot, but its length and amino acid composition are crucial for modulating the kinetics of channel inactivation.'
    }

    # Print the answer to the first question
    print("1) Hotspots on the beta-1 subunit for interaction with the alpha-1 subunit:\n" + "-"*80)
    print("The primary interaction is mediated by the Beta Interaction Domain (BID) within the GK domain of the beta-1 subunit. Key residues forming the binding pocket are:\n")
    for residue, description in interaction_hotspots.items():
        res_code, res_pos_type = residue.split(' ')
        res_pos = res_pos_type.strip('()')
        print(f"Residue: {res_code} at Position {res_pos} -> Role: {description}")

    print("\n\n" + "="*80 + "\n")

    # Print the answer to the second question
    print("2) Hotspots on the beta-1 subunit for fine-tuning alpha-1 subunit gating properties:\n" + "-"*80)
    print("Gating modulation is complex, involving regions outside the primary binding interface. Key residues and regions implicated in this regulation are:\n")
    for residue, description in gating_modulation_hotspots.items():
        # Handle the region entry differently from single residues
        if "Region" in residue:
            print(f"Region: {residue} -> Role: {description}")
        else:
            res_code, res_pos_type = residue.split(' ')
            res_pos = res_pos_type.strip('()')
            print(f"Residue: {res_code} at Position {res_pos} -> Role: {description}")


if __name__ == '__main__':
    get_calcium_channel_hotspots()
    # In a real biological context, you would now use this information for experiments,
    # such as site-directed mutagenesis to validate these functional roles.
    # The final deliverable is the printed output above.
    # The format required is just the list of amino acids and positions. Let's make sure the printout is very clear about that.

    print("\n\n<<<")
    print("Interaction Hotspots: Y253, M275, V278, W282, L312, A315, L319, I320, V323")
    print("Gating Modulation Hotspots: E134, E136, W161, and the HOOK Region (residues ~219-245)")
    print(">>>")
