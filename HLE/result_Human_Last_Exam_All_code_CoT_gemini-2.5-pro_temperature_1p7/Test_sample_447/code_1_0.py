import collections

def analyze_mouse_mutations():
    """
    Analyzes the effect of specific genetic mutations on the antibody response in mice.
    """
    # Define the groups and their associated genes/mutations
    groups = {
        'G1': 'AID-(V18R)',
        'G2': 'CD40-KO',
        'G3': 'H2-IAd-(E137A/V142A)',
        'G4': 'CD8-(V247D)',
        'G5': 'H2-IAd-(T139A)',
        'G6': 'MyD88-KO'
    }

    # Store explanations for why a group's antibody titer would be different
    explanations = {
        'G1': "AID is the essential enzyme for Somatic Hypermutation (SHM). A mutation in AID directly impairs the generation of hypermutated antibodies.",
        'G2': "CD40 on B cells must bind CD40L on T cells for T-cell help, which is required for germinal center formation and SHM. A knockout (KO) ablates this process.",
        'G3': "H2-IAd is the MHC-II molecule. Mutations here can impair antigen presentation to T helper cells, reducing the T-cell help needed for the antibody response.",
        'G4': "CD8 is on cytotoxic T cells, which are not directly involved in providing help to B cells for antibody production. This mutation is unlikely to affect the measured response.",
        'G5': "H2-IAd is the MHC-II molecule. Similar to G3, mutations here can impair antigen presentation to T helper cells, reducing T-cell help.",
        'G6': "MyD88 is the key adaptor for the CpG adjuvant (a TLR9 agonist). Its knockout prevents the adjuvant from boosting the immune response, leading to a weaker outcome than in wild-type."
    }

    print("Analysis of mutant mouse groups for expected differences in high-affinity antibody titers:")
    print("-" * 80)

    affected_groups = []
    for group_id, mutation in groups.items():
        if group_id != 'G4':
            affected_groups.append(group_id)
            print(f"Group {group_id} [{mutation}]: EXPECTED to be significantly different.")
            print(f"  Reason: {explanations[group_id]}\n")
        else:
            print(f"Group {group_id} [{mutation}]: NOT expected to be significantly different.")
            print(f"  Reason: {explanations[group_id]}\n")
            
    print("-" * 80)
    print("Conclusion: The groups in which the titer of high-affinity OVA-specific antibodies")
    print("that have undergone SHM would be expected to be significantly different are:")
    # Printing the numbers in the final selection
    final_group_numbers = [s.replace('G', '') for s in affected_groups]
    print(f"G{', G'.join(final_group_numbers)}")

analyze_mouse_mutations()
<<<C>>>