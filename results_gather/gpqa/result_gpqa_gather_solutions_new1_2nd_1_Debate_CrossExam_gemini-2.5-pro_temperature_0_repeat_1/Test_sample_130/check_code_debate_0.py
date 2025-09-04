def check_diels_alder_noesy():
    """
    This function checks the correctness of the answer to the Diels-Alder NOESY problem
    by codifying the chemical reasoning into a logical sequence.
    """

    # 1. Define the problem's chemical knowledge base
    # This includes proton assignments, stereochemical rules, and spatial proximities.
    chemical_knowledge = {
        "proton_assignments": {
            "H_anhydride": {"signal": "2H singlet", "shift_ppm": 3.5},
            "Me_vinyl": {"signal": "6H singlet", "shift_ppm": 1.7},
            "H_bridge": {"signal": "1H doublet", "shift_ppm": 1.5},
            "Me_bridgehead": {"signal": "6H singlet", "shift_ppm": 1.0}
        },
        "diene_properties": {
            "1,2,3,4-tetramethyl-1,3-cyclopentadiene": {
                "is_sterically_hindered": True
            }
        },
        "diels_alder_rules": {
            "default_major_product": "endo",  # The Alder-Endo Rule
            "steric_hindrance_major_product": "exo" # The exception
        },
        "noesy_proximities": {
            # Key distinguishing interaction in each isomer
            "endo": sorted(["H_anhydride", "H_bridge"]),
            "exo": sorted(["H_anhydride", "Me_vinyl"])
        }
    }

    # Map the answer choices to the types of interacting protons
    answer_choices = {
        "A": sorted(["H_anhydride", "H_bridge"]),
        "B": sorted(["H_anhydride", "Me_vinyl"]),
        "C": sorted(["Me_bridgehead", "H_bridge"]),
        "D": sorted(["Me_bridgehead", "Me_vinyl"])
    }

    # The final answer provided by the LLM to be checked
    provided_answer = "B"

    # 2. Execute the logical deduction process

    # Step A: Determine the major product based on steric hindrance
    diene = "1,2,3,4-tetramethyl-1,3-cyclopentadiene"
    if chemical_knowledge["diene_properties"][diene]["is_sterically_hindered"]:
        major_product_type = chemical_knowledge["diels_alder_rules"]["steric_hindrance_major_product"]
    else:
        major_product_type = chemical_knowledge["diels_alder_rules"]["default_major_product"]

    # Step B: Identify the unique NOESY interaction for the major product
    # This interaction is present in the major product but absent/weak in the minor.
    key_interaction = chemical_knowledge["noesy_proximities"][major_product_type]

    # Step C: Find which answer choice corresponds to this key interaction
    derived_correct_answer = None
    for option, interaction in answer_choices.items():
        if interaction == key_interaction:
            derived_correct_answer = option
            break

    if derived_correct_answer is None:
        return "Logic Error: The script could not find an answer choice matching the derived interaction."

    # 3. Compare the derived answer with the provided answer and return the result
    if derived_correct_answer == provided_answer:
        return "Correct"
    else:
        reason = (
            f"The provided answer '{provided_answer}' is incorrect. The logically derived answer is '{derived_correct_answer}'.\n\n"
            f"Reasoning:\n"
            f"1. Constraint Check (Stereochemistry): The diene (1,2,3,4-tetramethyl-1,3-cyclopentadiene) is exceptionally bulky. This steric hindrance overrides the standard Alder-Endo rule.\n"
            f"2. Conclusion (Major Product): The major product is the 'exo' adduct, not the 'endo' adduct.\n"
            f"3. Constraint Check (NOESY): The question asks for a NOESY cross-peak present in the major product. We must therefore analyze the 'exo' structure.\n"
            f"4. Analysis (Spatial Proximity): In the 'exo' adduct, the anhydride protons (~3.5 ppm) are spatially close to the vinylic methyl protons (~1.7 ppm).\n"
            f"5. Final Answer Derivation: This interaction corresponds to option '{derived_correct_answer}'.\n\n"
            f"The provided answer '{provided_answer}' is wrong because it either incorrectly assumes the 'endo' product is major or misidentifies the key spatial interaction."
        )
        return reason

# Run the check
result = check_diels_alder_noesy()
print(result)