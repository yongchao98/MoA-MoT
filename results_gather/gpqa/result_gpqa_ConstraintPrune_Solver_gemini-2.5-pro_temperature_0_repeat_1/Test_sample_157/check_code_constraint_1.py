def check_mutation_phenotype():
    """
    This function models the biological scenario from the question to logically
    verify the correctness of the provided answer.

    It is based on the following premises derived from the question:
    1.  Haplosufficiency: A cell with 50% of the normal amount of functional protein
        exhibits a wild-type (normal) phenotype. This is established by the
        recessive nature of mutation X.
    2.  Dominant-Negative Effect: Mutation Y is dominant-negative, meaning a
        heterozygous cell (WT/Y) has a loss-of-function phenotype.
    3.  Mechanism of Interference: For a dominant-negative effect to occur under
        haplosufficiency, the mutant protein (Y) must actively interfere with and
        reduce the function of the wild-type (WT) protein from the other allele.
    4.  Location: Mutation Y is in the dimerization domain, suggesting the
        interference happens at the level of protein-protein interaction.
    """

    # Let's define a function to check a proposed mechanism
    def evaluate_mechanism(option, description, expected_phenotype, mechanism_logic):
        print(f"--- Evaluating Option {option}: {description} ---")
        
        # All dominant-negative mutations are a class of loss-of-function.
        if "gain-of-function" in expected_phenotype:
            reason = "This is incorrect. A dominant-negative mutation is a type of loss-of-function, not gain-of-function."
            print(f"Result: Fails. {reason}\n")
            return False, reason

        # Simulate the mechanism
        final_phenotype = mechanism_logic()

        # Check if the resulting phenotype matches the dominant-negative requirement
        if final_phenotype == "loss-of-function":
            # The mechanism correctly produces a loss-of-function phenotype
            if "loss-of-function" in expected_phenotype:
                reason = "This mechanism is consistent with a dominant-negative effect under haplosufficiency."
                print(f"Result: Plausible. {reason}\n")
                return True, reason
            else:
                # This case is for options that predict a different phenotype but the mechanism leads to LoF
                reason = f"The mechanism leads to a loss-of-function phenotype, which contradicts the option's claim of a '{expected_phenotype}'."
                print(f"Result: Fails. {reason}\n")
                return False, reason
        elif final_phenotype == "wild-type":
            # The mechanism fails to produce a loss-of-function phenotype
            reason = "This mechanism results in a wild-type phenotype due to haplosufficiency, which contradicts the dominant-negative effect described in the question."
            print(f"Result: Fails. {reason}\n")
            return False, reason
        
        return False, "Unknown logic error."


    # --- Define the logic for each mechanism ---

    # In a WT/Y cell, 50% of protein is WT. Due to haplosufficiency, if these WT proteins
    # can form functional dimers freely, the phenotype will be wild-type.
    
    def mechanism_A_logic():
        # Y protein cannot dimerize. WT proteins are free to form WT-WT dimers.
        # This results in 50% of normal functional protein level.
        # Due to haplosufficiency, this is a wild-type phenotype.
        return "wild-type"

    def mechanism_C_logic():
        # Y protein is degraded. Only WT proteins remain.
        # This results in 50% of normal functional protein level.
        # Due to haplosufficiency, this is a wild-type phenotype.
        return "wild-type"

    def mechanism_D_logic():
        # Y protein dimerizes with WT protein, sequestering it into non-functional aggregates.
        # This reduces the amount of functional WT-WT dimers to well below 50% of normal.
        # For example, with random dimerization, only 25% of dimers are WT-WT.
        # This is insufficient for normal function, leading to a loss-of-function phenotype.
        return "loss-of-function"

    # --- Run the evaluation ---
    
    options = {
        'A': ("loss of protein dimerization and wild-type phenotype", mechanism_A_logic),
        'B': ("change of protein conformation and gain-of-function phenotype", lambda: "gain-of-function"),
        'C': ("protein degradation and loss-of-function of the wild-type allele", mechanism_C_logic),
        'D': ("protein aggregation and loss-of-function phenotype", mechanism_D_logic)
    }
    
    correct_options = []
    error_reasons = {}

    for option, (description, logic) in options.items():
        is_correct, reason = evaluate_mechanism(option, description, description, logic)
        if is_correct:
            correct_options.append(option)
        else:
            error_reasons[option] = reason

    # --- Final Verdict ---
    llm_answer = 'D'
    if llm_answer in correct_options and len(correct_options) == 1:
        return "Correct"
    elif llm_answer not in correct_options:
        return f"Incorrect. The answer 'D' is wrong because: {error_reasons.get(llm_answer, 'It fails the logical check.')}"
    else:
        return f"Incorrect. The answer 'D' is plausible, but other options ({correct_options}) are also plausible according to the model."


# Execute the check and print the result
result = check_mutation_phenotype()
print(f"\n-------------------\nFinal Check Result: {result}\n-------------------")
if result != "Correct":
    print(f"The provided answer is incorrect. Reason: {result}")
