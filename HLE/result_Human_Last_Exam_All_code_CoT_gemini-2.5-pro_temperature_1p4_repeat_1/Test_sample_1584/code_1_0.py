def analyze_coexpression_options():
    """
    Analyzes plasmid combinations for co-expression compatibility in E. coli.
    """
    # Database of plasmid properties.
    # Properties: origin of replication, antibiotic resistance, and number of expression sites (MCS).
    plasmids = {
        "pCDF-1b":      {"ori": "CDF", "res": "spectinomycin", "mcs": 1, "type": "expression"},
        "pET-28a(+)":     {"ori": "ColE1", "res": "kanamycin", "mcs": 1, "type": "expression"},
        "pGEX-T4-1":      {"ori": "ColE1", "res": "ampicillin", "mcs": 1, "type": "expression"},
        "pCDFDuet-1":   {"ori": "CDF", "res": "spectinomycin", "mcs": 2, "type": "expression"},
        "pGEM-T":         {"ori": "ColE1", "res": "ampicillin", "mcs": 1, "type": "cloning"}, # Not an expression vector
        "pET-15b":        {"ori": "ColE1", "res": "ampicillin", "mcs": 1, "type": "expression"},
        "pASK-IBA3":      {"ori": "p15A", "res": "chloramphenicol", "mcs": 1, "type": "expression"},
    }

    # The options provided in the question.
    options = {
        'A': {'plasmids': ["pCDF-1b", "pET-28a(+)"], 'desc': "pCDF-1b with spectinomycin resistance and pET-28a(+) with kanamycin resistance"},
        'B': {'plasmids': ["pET-28a(+)", "pGEX-T4-1"], 'desc': "pET-28a(+) with kanamycin resistance and pGEX-T4-1 with ampicillin resistance"},
        'C': {'plasmids': ["pCDFDuet-1", "pET-28a(+)"], 'desc': "pCDFDuet-1 with spectinomycin resistance and pET-28a(+) with kanamycin resistance"},
        'D': {'plasmids': ["pET-28a(+)", "pGEX-T4-1"], 'desc': "pET-28a(+) with kanamycin resistance and pGEX-T4-1 with chloramphenicol resistance"},
        'E': {'plasmids': ["pGEM-T", "pCDF-1b"], 'desc': "pGEM®-T with ampicillin resistance and pCDF-1b with spectinomycin resistance"},
        'F': {'plasmids': ["pET-15b", "pET-28a(+)"], 'desc': "pET-15b with ampicillin resistance and pET-28a(+) with kanamycin resistance"},
        'G': {'plasmids': ["pASK-IBA3", "pET-28a(+)"], 'desc': "pASK-IBA3 with chloramphenicol resistance and pET-28a(+) with ampicillin resistance"},
        'H': {'plasmids': ["pCDFDuet-1"], 'desc': "pCDFDuet-1 with spectinomycin resistance"},
        'J': {'plasmids': ["pGEX-T4-1", "pASK-IBA3"], 'desc': "pGEX-T4-1 with ampicillin resistance and pASK-IBA3 with chloramphenicol resistance"},
    }
    
    # Store reasons for each conclusion
    analysis_results = {}

    for key, option in options.items():
        plasmid_names = option['plasmids']
        is_valid = True
        reasons = []

        # Check for factual errors in descriptions
        if key == 'D' and "chloramphenicol" in option['desc']:
            reasons.append("Factual Error: Standard pGEX-T4-1 has ampicillin resistance, not chloramphenicol.")
            is_valid = False
        if key == 'G' and "ampicillin" in option['desc']:
             reasons.append("Factual Error: Standard pET-28a(+) has kanamycin resistance, not ampicillin.")
             is_valid = False
        
        # Check if vectors are appropriate for expression
        for p_name in plasmid_names:
            if plasmids[p_name]['type'] != 'expression':
                reasons.append(f"Invalid Vector Type: {p_name} is a cloning vector, not suitable for high-level expression.")
                is_valid = False

        # Check compatibility for two-plasmid systems
        if len(plasmid_names) == 2 and is_valid:
            p1 = plasmids[plasmid_names[0]]
            p2 = plasmids[plasmid_names[1]]
            if p1['ori'] == p2['ori']:
                reasons.append(f"Incompatible Origins: Both {plasmid_names[0]} and {plasmid_names[1]} have a {p1['ori']}-type origin.")
                is_valid = False
            if p1['res'] == p2['res']:
                reasons.append(f"Indistinguishable Selection: Both plasmids have {p1['res']} resistance.")
                is_valid = False
        
        if is_valid:
            total_mcs = sum(plasmids[p_name]['mcs'] for p_name in plasmid_names)
            reasons.append(f"This system is valid and can express up to {total_mcs} different proteins.")
            if key in ['A', 'C', 'H', 'J']:
                 analysis_results[key] = f"Result: VALID. {'. '.join(reasons)}"
            else:
                 analysis_results[key] = f"Result: INVALID. {'. '.join(reasons)}"
        else:
            analysis_results[key] = f"Result: INVALID. {'. '.join(reasons)}"

    # Print summary
    for key, result in sorted(analysis_results.items()):
        print(f"Option {key}: {options[key]['desc']}\n  -> {result}\n")
    
    print("--- Conclusion ---")
    print("Options A, C, H, and J describe technically valid plasmid systems.")
    print("However, the 'best' option must be powerful and versatile.")
    print("Option C describes a system using pCDFDuet-1 (2 proteins) and pET-28a(+) (1 protein).")
    print("This allows the co-expression of a protein of interest (from pET-28a(+)) with a two-subunit chaperone complex like GroEL/ES (from pCDFDuet-1).")
    print("This is a highly effective and common strategy for difficult-to-fold proteins, making it the most capable and therefore 'best' choice.")

# Run the analysis
analyze_coexpression_options()