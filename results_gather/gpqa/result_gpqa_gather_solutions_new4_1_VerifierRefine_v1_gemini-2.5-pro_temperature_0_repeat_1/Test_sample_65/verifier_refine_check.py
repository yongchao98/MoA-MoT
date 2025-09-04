def check_smeft_answer(llm_answer: str):
    """
    Checks the correctness of an answer to the SMEFT symmetries question.

    The question asks which symmetries must be respected by all operators in the SMEFT.
    1. Lorentz Symmetry
    2. Poincare symmetry
    3. CP symmetry
    4. CPT symmetry

    Args:
        llm_answer: The letter ('A', 'B', 'C', or 'D') corresponding to the chosen answer.

    Returns:
        A string indicating "Correct" or the reason for the incorrectness.
    """
    # Define the ground truth based on established principles of Quantum Field Theory
    # SMEFT is a local, relativistic QFT, so it must respect Poincare symmetry (which includes Lorentz)
    # and CPT symmetry (by the CPT theorem).
    # CP is violated in the SM and SMEFT is designed to parameterize new sources of CP violation,
    # so it is NOT a required symmetry.
    correct_symmetries = {1, 2, 4}
    
    symmetry_names = {
        1: "Lorentz Symmetry",
        2: "Poincare symmetry",
        3: "CP symmetry",
        4: "CPT symmetry"
    }

    # Define the options as provided in the question
    options = {
        'A': {1, 2},
        'B': {3, 4},
        'C': {1, 2, 4},
        'D': {1, 3, 4}
    }

    # Sanitize input
    answer_key = llm_answer.strip().upper()

    if answer_key not in options:
        return f"Invalid answer format. The provided answer was '{llm_answer}', but it must be one of 'A', 'B', 'C', or 'D'."

    proposed_symmetries = options[answer_key]

    # Check if the proposed set of symmetries matches the correct set
    if proposed_symmetries == correct_symmetries:
        return "Correct"
    else:
        # Analyze the discrepancy
        missing_symmetries = correct_symmetries - proposed_symmetries
        incorrectly_included_symmetries = proposed_symmetries - correct_symmetries
        
        reasons = []
        if incorrectly_included_symmetries:
            names = ", ".join([symmetry_names[s] for s in sorted(list(incorrectly_included_symmetries))])
            reasons.append(f"The answer incorrectly includes symmetries that are not required for all SMEFT operators: {names}. CP symmetry is not a required symmetry.")
        
        if missing_symmetries:
            names = ", ".join([symmetry_names[s] for s in sorted(list(missing_symmetries))])
            reasons.append(f"The answer is missing required symmetries: {names}. All SMEFT operators must respect Lorentz, Poincare, and CPT symmetry.")
            
        return "Incorrect. " + " ".join(reasons)

# The final answer provided by the LLM was <<<C>>>.
# Let's check it.
final_answer_from_llm = "C"
result = check_smeft_answer(final_answer_from_llm)
print(result)

# You can also test other options to see the failure reasons.
# print("\n--- Testing other options ---")
# print(f"Checking A: {check_smeft_answer('A')}")
# print(f"Checking B: {check_smeft_answer('B')}")
# print(f"Checking D: {check_smeft_answer('D')}")
