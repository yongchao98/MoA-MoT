def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer by modeling the problem's constraints.

    The function analyzes the biological principles described in the question:
    1.  Haplosufficiency: Established by the recessive nature of mutation X, meaning one
        functional allele is sufficient for a wild-type phenotype.
    2.  Dominant-Negative Effect: Mutation Y causes a loss-of-function phenotype in a
        heterozygote by interfering with the wild-type protein.
    3.  Mechanism of Interference: The mutation is in the dimerization domain, so the
        interference must be related to protein-protein interaction.
    """

    # --- Define Constraints from the Question ---

    # A dominant-negative mutation results in a loss-of-function phenotype.
    EXPECTED_PHENOTYPE = "loss_of_function"

    # The core principle of a dominant-negative mutation is that the mutant protein
    # interferes with the function of the wild-type protein.
    MUST_INTERFERE_WITH_WT = True

    # To interfere, the mutant protein must be present and stable. If it's degraded,
    # it cannot interact with the wild-type protein.
    MUTANT_PROTEIN_MUST_BE_PRESENT = True
    
    # The mutation is in the dimerization domain, so the proposed molecular mechanism
    # must be a plausible consequence of a faulty protein-protein interaction domain.
    MECHANISM_MUST_BE_PLAUSIBLE = True

    # --- Define Properties of Each Option ---
    options = {
        "A": {
            "phenotype": "wild_type",
            "interferes": False,  # If Y can't dimerize, it can't interfere with WT via dimerization.
            "mutant_present": True,
            "mechanism_plausible": True
        },
        "B": {
            "phenotype": "gain_of_function",
            "interferes": True,
            "mutant_present": True,
            "mechanism_plausible": True
        },
        "C": {
            "phenotype": "loss_of_function",
            "interferes": False,  # A degraded protein cannot interfere.
            "mutant_present": False,
            "mechanism_plausible": False # Degradation is a consequence, not a direct interference mechanism.
        },
        "D": {
            "phenotype": "loss_of_function",
            "interferes": True,   # Aggregation traps WT protein, which is a form of interference.
            "mutant_present": True,
            "mechanism_plausible": True # Aggregation is a plausible result of a mutation in a dimerization domain.
        }
    }

    llm_answer_choice = "D"
    selected_option = options[llm_answer_choice]

    # --- Perform Checks ---

    # 1. Check Phenotype Constraint
    if selected_option["phenotype"] != EXPECTED_PHENOTYPE:
        return (f"Incorrect. The question describes a dominant-negative mutation, which by definition "
                f"causes a '{EXPECTED_PHENOTYPE}' phenotype. The answer proposes a '{selected_option['phenotype']}' phenotype.")

    # 2. Check Mutant Protein Presence Constraint
    if selected_option["mutant_present"] != MUTANT_PROTEIN_MUST_BE_PRESENT:
        return (f"Incorrect. The proposed mechanism (protein degradation) means the mutant protein is not present to "
                f"interfere with the wild-type protein. Without the mutant protein, the cell would have a wild-type "
                f"phenotype due to haplosufficiency, which contradicts the dominant-negative effect.")

    # 3. Check Interference Constraint
    if selected_option["interferes"] != MUST_INTERFERE_WITH_WT:
        return (f"Incorrect. A dominant-negative mutation requires the mutant protein to interfere with the wild-type "
                f"protein. The proposed mechanism does not involve interference. This would lead to a wild-type "
                f"phenotype due to haplosufficiency, not the observed loss-of-function.")

    # 4. Check Mechanism Plausibility Constraint
    if selected_option["mechanism_plausible"] != MECHANISM_MUST_BE_PLAUSIBLE:
        return (f"Incorrect. The proposed mechanism is not a plausible direct consequence of a mutation in the "
                f"dimerization domain that leads to interference. Protein degradation is a cellular response, not "
                f"a mechanism of interference itself.")

    # If all constraints are satisfied, the answer is correct.
    return "Correct"

# Run the check
result = check_answer_correctness()
print(result)