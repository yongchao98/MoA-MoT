def check_molecular_phenotype_answer(final_answer):
    """
    Checks the correctness of the answer for the dominant-negative mutation question.

    The function codifies the following biological principles:
    1.  A dominant-negative mutation is a type of loss-of-function, not gain-of-function.
    2.  A dominant mutation causes a mutant phenotype, not a wild-type one.
    3.  For a dominant-negative effect in a dimeric protein, the mutant must be able to
        interact (dimerize) with the wild-type protein to interfere. A complete loss
        of dimerization would lead to a recessive phenotype.
    4.  The most accurate description of the dominant-negative effect is the functional
        inactivation of the wild-type allele's product.
    5.  A common cellular mechanism for handling faulty protein complexes (like a
        poisoned dimer) is targeting them for degradation.
    """

    # Define the properties of each option
    options = {
        'A': {
            "mechanism": "loss of protein dimerization",
            "phenotype": "wild-type phenotype",
            "is_correct": False,
            "reason": "This is incorrect. A dominant mutation causes a mutant phenotype, not wild-type. Also, if the mutant protein could not dimerize, it could not interfere with the wild-type protein, making the mutation recessive, not dominant-negative."
        },
        'B': {
            "mechanism": "protein aggregation",
            "phenotype": "loss-of-function phenotype",
            "is_correct": False,
            "reason": "This is plausible but less precise than option D. While 'loss-of-function' is correct, 'protein aggregation' is a general mechanism. Option D provides a more specific and common mechanism (degradation) and a more accurate description of the genetic effect ('loss-of-function of the wild-type allele')."
        },
        'C': {
            "mechanism": "change of protein conformation",
            "phenotype": "gain-of-function phenotype",
            "is_correct": False,
            "reason": "This is incorrect. A dominant-negative mutation causes a loss of function, not a gain of function."
        },
        'D': {
            "mechanism": "protein degradation",
            "phenotype": "loss-of-function of the wild-type allele",
            "is_correct": True,
            "reason": "This is the most accurate answer. 'Loss-of-function of the wild-type allele' is the definition of a dominant-negative effect. 'Protein degradation' is a common and plausible molecular mechanism for how a cell handles the resulting non-functional heterodimers."
        }
    }

    # Check the provided final_answer
    if final_answer not in options:
        return f"Invalid answer choice '{final_answer}'. Please provide one of {list(options.keys())}."

    selected_option = options[final_answer]

    if selected_option["is_correct"]:
        return "Correct"
    else:
        return f"Incorrect. The answer '{final_answer}' is wrong because: {selected_option['reason']}"

# The final answer provided by the LLM was 'D'.
llm_final_answer = 'D'

# Run the check
result = check_molecular_phenotype_answer(llm_final_answer)
print(result)