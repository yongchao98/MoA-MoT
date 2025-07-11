def find_affected_groups():
    """
    Analyzes the function of genes mutated in different mouse groups
    to determine which would have an altered high-affinity antibody response.
    """
    
    # Description of each group's gene function in the context of antibody response
    group_analysis = {
        'G1 (AID)': "AID is essential for somatic hypermutation (SHM). Mutation will impair affinity maturation.",
        'G2 (CD40)': "CD40 is crucial for T-cell help and germinal center formation. KO will prevent affinity maturation.",
        'G3 (H2-IAd)': "MHC Class II is essential for presenting antigen to helper T-cells. Mutation can block T-cell help.",
        'G4 (CD8)': "CD8 is on cytotoxic T-cells, not directly involved in helper T-cell-mediated B-cell help for antibody production.",
        'G5 (H2-IAd)': "MHC Class II is essential for presenting antigen to helper T-cells. Mutation can block T-cell help.",
        'G6 (MyD88)': "MyD88 is essential for the CpG adjuvant effect via TLR9. KO will result in a much weaker immune response."
    }

    # Criteria for significant difference: The gene must be crucial for
    # T-dependent B-cell affinity maturation or for the adjuvant effect.
    affected_groups = []
    
    # CD8 (G4) is the only gene not directly involved in this pathway.
    unaffected_gene_symbol = "CD8"
    
    for group, description in group_analysis.items():
        if unaffected_gene_symbol not in group:
            # Extract just the group number, e.g., 'G1'
            affected_groups.append(group.split(' ')[0])

    print("The groups expected to have a significantly different titer of high-affinity OVA-specific antibodies are:")
    # We join the list for a clean output.
    # The final answer format is G1, G2, G3, G5, G6 which corresponds to choice C.
    print(', '.join(affected_groups))

find_affected_groups()