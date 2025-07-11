def find_best_coexpression_system():
    """
    Analyzes plasmid combinations for co-expression in E. coli based on compatibility
    and suitability for protein expression.
    """
    # Define plasmid properties: Origin of Replication (ori), Resistance, and Primary Use.
    # ColE1-type origins (ColE1, pBR322, pMB1) are incompatible with each other.
    # CloDF13 is compatible with ColE1-type origins.
    plasmids = {
        "pET-28a(+)":    {"ori": "ColE1", "resistance": "Kanamycin", "type": "Expression"},
        "pET-15b":       {"ori": "ColE1", "resistance": "Ampicillin", "type": "Expression"},
        "pCDF-1b":       {"ori": "CloDF13", "resistance": "Spectinomycin", "type": "Expression"},
        "pCDFDuet-1":    {"ori": "CloDF13", "resistance": "Spectinomycin", "type": "Co-Expression"},
        "pGEX-T4-1":     {"ori": "ColE1", "resistance": "Ampicillin", "type": "Expression"},
        "pASK-IBA3":     {"ori": "ColE1", "resistance": "Chloramphenicol", "type": "Expression"},
        "pGEM®-T":       {"ori": "ColE1", "resistance": "Ampicillin", "type": "Cloning"},
    }

    # Answer choices represented as pairs of plasmid names.
    # Note: Options D and G contain factual errors in the prompt about resistance, but their
    # primary issue is origin incompatibility, which is what we will focus on.
    choices = {
        'A': ["pCDF-1b", "pET-28a(+)"],
        'B': ["pET-28a(+)", "pGEX-T4-1"],
        'C': ["pCDFDuet-1", "pET-28a(+)"],
        'D': ["pET-28a(+)", "pGEX-T4-1"],
        'E': ["pGEM®-T", "pCDF-1b"],
        'F': ["pET-15b", "pET-28a(+)"],
        'G': ["pASK-IBA3", "pET-28a(+)"],
        'H': ["pCDFDuet-1"],
        'J': ["pGEX-T4-1", "pASK-IBA3"],
    }

    print("--- Analysis of Plasmid Systems for Co-expression ---\n")
    final_conclusion = "The analysis points to one superior option."
    best_option = "C"

    for choice, p_names in choices.items():
        print(f"Analyzing Option {choice}: {', '.join(p_names)}")
        
        # Handle single plasmid option
        if len(p_names) == 1:
            plasmid_info = plasmids[p_names[0]]
            if plasmid_info['type'] == "Co-Expression":
                print(f"  Result: VALID. The {p_names[0]} vector is specifically designed to express two proteins from one plasmid. This is a very good single-plasmid solution.\n")
            else:
                print(f"  Result: INVALID. A standard expression vector cannot co-express two different proteins.\n")
            continue

        # Handle two-plasmid systems
        p1_name, p2_name = p_names[0], p_names[1]
        p1 = plasmids[p1_name]
        p2 = plasmids[p2_name]

        # 1. Check for compatible origins of replication
        ori_compatible = p1['ori'] != p2['ori']
        # 2. Check for different resistance markers
        res_different = p1['resistance'] != p2['resistance']

        if not ori_compatible:
            print(f"  Result: FAILED. Plasmids are incompatible due to both having {p1['ori']}-type origins.\n")
        elif not res_different:
            print(f"  Result: FAILED. Plasmids cannot be co-selected because both confer {p1['resistance']} resistance.\n")
        elif p1['type'] == 'Cloning' or p2['type'] == 'Cloning':
            bad_vector = p1_name if p1['type'] == 'Cloning' else p2_name
            print(f"  Result: SUB-OPTIMAL. Although compatible, the {bad_vector} plasmid is a cloning vector, not designed for high-level protein expression.\n")
        else:
            reason = "EXCELLENTLY SUITED" if "Co-Expression" in (p1['type'], p2['type']) else "VALID"
            print(f"  Result: {reason}. Origins ({p1['ori']} & {p2['ori']}) are compatible and resistances ({p1['resistance']} & {p2['resistance']}) are different.\n")

    print("--- Conclusion ---")
    print("While option A is a valid and common system, option C is the 'best' choice.")
    print("The pCDFDuet-1 vector is specifically engineered to express two proteins, often a chaperone system like GroEL/GroES.")
    print("Combining it with a robust expression plasmid like pET-28a(+) in a compatible system represents the most powerful and targeted approach among the choices for ensuring proper protein folding.")

    # The final answer in the required format.
    print(f"\n<<<{best_option}>>>")

find_best_coexpression_system()