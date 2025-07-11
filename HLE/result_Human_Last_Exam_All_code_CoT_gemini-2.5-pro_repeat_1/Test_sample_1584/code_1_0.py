def solve_coexpression_problem():
    """
    Analyzes plasmid combinations for co-expression in E. coli to find the best option.
    """
    # Database of plasmid properties
    # Properties: {origin group, standard resistance, vector type}
    plasmids_db = {
        "pCDF-1b": {"origin": "CloDF13", "resistance": "spectinomycin", "type": "expression"},
        "pET-28a(+)": {"origin": "ColE1", "resistance": "kanamycin", "type": "expression"},
        "pGEX-T4-1": {"origin": "ColE1", "resistance": "ampicillin", "type": "expression"},
        "pCDFDuet-1": {"origin": "CloDF13", "resistance": "spectinomycin", "type": "co-expression"},
        "pGEM-T": {"origin": "ColE1", "resistance": "ampicillin", "type": "cloning"},
        "pET-15b": {"origin": "ColE1", "resistance": "ampicillin", "type": "expression"},
        "pASK-IBA3": {"origin": "ColE1", "resistance": "chloramphenicol", "type": "expression"},
    }

    # Options as defined in the question
    # Format: {Option_Letter: [(plasmid_name, resistance_in_question), ...]}
    options = {
        "A": [("pCDF-1b", "spectinomycin"), ("pET-28a(+)", "kanamycin")],
        "B": [("pET-28a(+)", "kanamycin"), ("pGEX-T4-1", "ampicillin")],
        "C": [("pCDFDuet-1", "spectinomycin"), ("pET-28a(+)", "kanamycin")],
        "D": [("pET-28a(+)", "kanamycin"), ("pGEX-T4-1", "chloramphenicol")],
        "E": [("pGEM-T", "ampicillin"), ("pCDF-1b", "spectinomycin")],
        "F": [("pET-15b", "ampicillin"), ("pET-28a(+)", "kanamycin")],
        "G": [("pASK-IBA3", "chloramphenicol"), ("pET-28a(+)", "ampicillin")],
        "H": [("pCDFDuet-1", "spectinomycin")],
        "J": [("pGEX-T4-1", "ampicillin"), ("pASK-IBA3", "chloramphenicol")],
    }

    print("Analyzing Plasmid Co-Expression Options:\n")
    
    valid_options = {}

    for letter, plasmids in sorted(options.items()):
        analysis = []
        is_valid = True
        
        if len(plasmids) == 2:
            p1_name, r1_q = plasmids[0]
            p2_name, r2_q = plasmids[1]
            p1_info = plasmids_db.get(p1_name)
            p2_info = plasmids_db.get(p2_name)

            if p1_info["origin"] == p2_info["origin"]:
                analysis.append(f"FAIL: Incompatible origins of replication (both are {p1_info['origin']})")
                is_valid = False
            else:
                analysis.append(f"PASS: Compatible origins ({p1_info['origin']} and {p2_info['origin']})")

            if r1_q == r2_q:
                analysis.append(f"FAIL: Identical selection markers ({r1_q})")
                is_valid = False
            else:
                analysis.append(f"PASS: Different selection markers ({r1_q} and {r2_q})")
            
            if p1_info["type"] == "cloning" or p2_info["type"] == "cloning":
                analysis.append("FAIL: Contains a cloning vector, not suitable for expression.")
                is_valid = False
            else:
                analysis.append("PASS: Both are expression-capable vectors.")

        elif len(plasmids) == 1:
            p_name, r_q = plasmids[0]
            p_info = plasmids_db.get(p_name)
            if p_info["type"] == "co-expression":
                analysis.append("PASS: A single co-expression vector is a valid strategy.")
                is_valid = True
            else:
                analysis.append("FAIL: Single plasmid is not a co-expression vector.")
                is_valid = False

        print(f"Option {letter}: {' '.join(p[0] for p in plasmids)}")
        for item in analysis:
            print(f"  - {item}")
        if is_valid:
            valid_options[letter] = True
            print("  --> Conclusion: VALID SYSTEM\n")
        else:
            print("  --> Conclusion: INVALID SYSTEM\n")

    print("--- Final Verdict ---")
    print("The valid systems are: " + ", ".join(valid_options.keys()))
    print("\nComparing the valid options (A, C, H):")
    print(" - Option A (pCDF-1b + pET-28a(+)): A classic and functional two-plasmid system.")
    print(" - Option H (pCDFDuet-1 alone): A robust single-plasmid system for co-expression.")
    print(" - Option C (pCDFDuet-1 + pET-28a(+)): The most powerful and versatile system.")
    print("\nOption C is the best choice because it combines a high-copy pET vector (for strong expression of the protein of interest) with a compatible pCDFDuet-1 vector. The 'Duet' vector is specifically designed for co-expression and can carry multiple chaperone genes (e.g., the GroEL/GroES complex), which is often required for effective protein folding assistance. This makes the system both powerful and adaptable.")

solve_coexpression_problem()
<<<C>>>