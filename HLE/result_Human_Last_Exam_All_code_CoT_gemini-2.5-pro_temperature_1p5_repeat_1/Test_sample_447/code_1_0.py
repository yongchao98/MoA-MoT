def analyze_mutant_groups():
    """
    Analyzes which mutant mouse groups would show a significantly different
    titer of high-affinity, somatically hypermutated (SHM) antibodies.
    """
    
    # Define the groups and the rationale for their effect on the antibody response
    groups_analysis = {
        'G1': ("AID-(V18R)", "AID is essential for SHM. A mutation here will impair affinity maturation.", True),
        'G2': ("CD40-KO", "CD40 is critical for T-cell help to B-cells and germinal center formation.", True),
        'G3': ("H2-IAd-(E137A/V142A)", "MHC-II mutation can impair antigen presentation to T-helper cells.", True),
        'G4': ("CD8-(V247D)", "CD8 is not central to T-helper cell-driven antibody production.", False),
        'G5': ("H2-IAd-(T139A)", "MHC-II mutation can impair antigen presentation to T-helper cells.", True),
        'G6': ("MyD88-KO", "MyD88 is required for the CpG adjuvant effect, which boosts the response.", True)
    }

    print("Analysis of each mutant group's expected antibody response:")
    
    affected_groups = []
    for group_id, (mutation, reason, is_affected) in groups_analysis.items():
        if is_affected:
            status = "AFFECTED"
            affected_groups.append(group_id)
        else:
            status = "NOT AFFECTED"
        
        print(f"- {group_id} [{mutation}]: {status}. Reason: {reason}")
        
    print("\n------------------------------------------------------")
    print("The groups expected to show a significantly different antibody titer are:")
    
    # The prompt asks to output each number/name in the final selection
    # We will format it as a list for clarity.
    final_selection_str = ", ".join(affected_groups)
    print(final_selection_str)

# Run the analysis
analyze_mutant_groups()
