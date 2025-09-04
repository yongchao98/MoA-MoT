def check_answer_correctness():
    """
    Checks the correctness of the final answer by verifying fundamental chemical constraints.
    """
    # --- Data Representation ---
    # Molecular formulas and functional groups of relevant compounds
    COMPOUNDS = {
        # Reactants
        "reactant_B": {"formula": "C8H10", "type": "Diyne"},
        "reactant_C": {"formula": "C7H12O", "type": "Allyl Vinyl Ether"},
        
        # Products for A
        "(Z)-1-(but-2-en-2-yloxy)ethen-1-amine": {"formula": "C6H11NO", "type": "Enamine Enol Ether"},
        "6-methyl-3,4-dihydro-2H-pyran-2-amine": {"formula": "C6H11NO", "type": "Cyclic Hemiaminal Ether"},

        # Products for B
        "(3Z,4E)-3,4-diethylidenecyclobut-1-ene": {"formula": "C8H10", "type": "Cyclobutene Derivative"},
        "(1Z,2E)-1,2-diethylidenecyclobutane": {"formula": "C8H12", "type": "Cyclobutane Derivative"},

        # Products for C
        "4-methylenehexanal": {"formula": "C7H12O", "type": "Carbonyl"},
        "4-methylenehexan-1-ol": {"formula": "C7H14O", "type": "Alcohol"},
    }

    # The products proposed by each multiple-choice option
    OPTIONS = {
        "A": {
            "A": "(Z)-1-(but-2-en-2-yloxy)ethen-1-amine",
            "B": "(3Z,4E)-3,4-diethylidenecyclobut-1-ene",
            "C": "4-methylenehexanal"
        },
        "B": {
            "A": "6-methyl-3,4-dihydro-2H-pyran-2-amine",
            "B": "(1Z,2E)-1,2-diethylidenecyclobutane",
            "C": "4-methylenehexanal"
        },
        "C": {
            "A": "(Z)-1-(but-2-en-2-yloxy)ethen-1-amine",
            "B": "(3Z,4E)-3,4-diethylidenecyclobut-1-ene",
            "C": "4-methylenehexan-1-ol"
        },
        "D": {
            "A": "6-methyl-3,4-dihydro-2H-pyran-2-amine",
            "B": "(1Z,2E)-1,2-diethylidenecyclobutane",
            "C": "4-methylenehexan-1-ol"
        }
    }

    llm_answer = "A"
    chosen_option = OPTIONS.get(llm_answer)

    if not chosen_option:
        return f"Invalid answer option '{llm_answer}'. Must be one of {list(OPTIONS.keys())}."

    errors = []

    # --- Constraint Check 1: Reaction B (Isomerization) ---
    reactant_b_formula = COMPOUNDS["reactant_B"]["formula"]
    product_b_name = chosen_option["B"]
    product_b_formula = COMPOUNDS[product_b_name]["formula"]

    if reactant_b_formula != product_b_formula:
        errors.append(
            f"Constraint failed for Reaction B: Product B '{product_b_name}' has formula {product_b_formula}, "
            f"which does not match the reactant's formula {reactant_b_formula}."
        )

    # --- Constraint Check 2: Reaction C (Claisen Rearrangement) ---
    reactant_c_formula = COMPOUNDS["reactant_C"]["formula"]
    product_c_name = chosen_option["C"]
    product_c_info = COMPOUNDS[product_c_name]
    product_c_formula = product_c_info["formula"]
    product_c_type = product_c_info["type"]

    # Check 2a: Isomerization
    if reactant_c_formula != product_c_formula:
        errors.append(
            f"Constraint failed for Reaction C: Product C '{product_c_name}' has formula {product_c_formula}, "
            f"which does not match the reactant's formula {reactant_c_formula}."
        )
    
    # Check 2b: Functional Group Type
    if product_c_type != "Carbonyl":
        errors.append(
            f"Constraint failed for Reaction C: The Claisen rearrangement must produce a carbonyl compound, "
            f"but the proposed product C '{product_c_name}' is an {product_c_type.lower()}."
        )

    # --- Final Verdict ---
    if not errors:
        return "Correct"
    else:
        return "\n".join(errors)

# Run the check
result = check_answer_correctness()
print(result)