def analyze_antibody_response():
    """
    Analyzes which mutant mouse strains would show a significantly different titer
    of high-affinity, somatically hypermutated OVA-specific antibodies compared to wild-type.
    """
    
    # A dictionary mapping each group to its mutated gene and role in immunity.
    group_info = {
        'G1': 'AID-(V18R): AID is the enzyme essential for Somatic Hypermutation (SHM).',
        'G2': 'CD40-KO: CD40 is critical for receiving T-cell help, essential for germinal center reactions.',
        'G3': 'H2-IAd-(E137A/V142A): H2-IAd (MHC-II) is required for antigen presentation to helper T-cells.',
        'G4': 'CD8-(V247D): CD8 is on cytotoxic T-cells, not directly involved in providing B-cell help.',
        'G5': 'H2-IAd-(T139A): H2-IAd (MHC-II) mutation, similar to G3, can impair T-cell help.',
        'G6': 'MyD88-KO: MyD88 is an essential adaptor for the CpG adjuvant (TLR9 agonist) signaling.'
    }

    # Groups expected to show a significant difference from wild-type.
    # We exclude G4 as CD8 T-cells are not directly involved in the B-cell help process for antibody production.
    affected_groups = ['G1', 'G2', 'G3', 'G5', 'G6']

    print("Analysis of Mutant Groups for Altered Antibody Response:")
    print("-" * 55)

    for group, reason in group_info.items():
        status = "AFFECTED" if group in affected_groups else "UNAFFECTED"
        print(f"Group {group} [{status}]: {reason}")
    
    print("-" * 55)
    print("The final list of groups expected to have a significantly different antibody titer are:")
    
    # Constructing the "equation" as requested, showing the numbers from the selected groups
    group_numbers = [group[1] for group in affected_groups]
    equation_str = " + ".join(f"G{num}" for num in group_numbers)
    print(f"Selected Groups = {equation_str}")


if __name__ == '__main__':
    analyze_antibody_response()
