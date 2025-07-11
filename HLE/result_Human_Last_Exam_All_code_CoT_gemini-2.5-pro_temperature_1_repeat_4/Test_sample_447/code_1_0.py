def analyze_antibody_response():
    """
    Analyzes which mutant mouse groups would have a significantly different
    titer of high-affinity, somatically hypermutated (SHM) antibodies.
    """
    mutant_groups = {
        'G1': {'gene': 'AID-(V18R)', 'role': 'Enzyme essential for Somatic Hypermutation (SHM).', 'affected': True},
        'G2': {'gene': 'CD40-KO', 'role': 'Receptor essential for T-cell help and Germinal Center (GC) formation.', 'affected': True},
        'G3': {'gene': 'H2-IAd-(E137A/V142A)', 'role': 'MHC Class II molecule, required for antigen presentation to helper T-cells.', 'affected': True},
        'G4': {'gene': 'CD8-(V247D)', 'role': 'Co-receptor on cytotoxic T-cells, not directly involved in helper T-cell/B-cell interaction for antibody production.', 'affected': False},
        'G5': {'gene': 'H2-IAd-(T139A)', 'role': 'MHC Class II molecule, required for antigen presentation to helper T-cells.', 'affected': True},
        'G6': {'gene': 'MyD88-KO', 'role': 'Adaptor for CpG adjuvant (TLR9) signaling, which enhances the immune response.', 'affected': True}
    }

    print("Analysis of Mutant Groups for High-Affinity Antibody Production:\n")
    
    affected_groups_numbers = []

    for group, data in mutant_groups.items():
        print(f"Group: {group} [{data['gene']}]")
        print(f"  - Role: {data['role']}")
        if data['affected']:
            print("  - Expected Outcome: Significantly different antibody titer compared to wild-type.")
            group_number = int(group.replace('G', ''))
            affected_groups_numbers.append(group_number)
        else:
            print("  - Expected Outcome: No significant difference in antibody titer compared to wild-type.")
        print("-" * 30)

    # Sort the numbers for a clean output
    affected_groups_numbers.sort()
    
    # The final "equation" is the list of affected group numbers.
    print("\nSummary:")
    print("The groups expected to show a significantly different titer of high-affinity, SHM antibodies are:")
    # Output each number in the final list
    final_equation = ", ".join(map(str, affected_groups_numbers))
    print(f"G{', G'.join(map(str, affected_groups_numbers))}")


if __name__ == '__main__':
    analyze_antibody_response()
    print("\n<<<C>>>")
