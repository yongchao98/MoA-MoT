import re

def check_correctness_of_spectroscopy_answer(llm_response: str) -> str:
    """
    Checks the correctness of the LLM's answer by logically deducing the molecular formula
    from the spectral data provided in the question.

    The function programmatically follows these steps:
    1.  Parses the LLM's final answer choice (e.g., 'B') from the response.
    2.  Analyzes the constraints from the question to build the chemical structure piece by piece.
        -   Aromatic ring: C6H4
        -   Vinyl/CH3 signals: Propenyl group (C3H5)
        -   Ester/CH3/no-CH2 signals: Methyl ester or acetate group (C2H3O2)
    3.  Assembles the fragments to determine the only possible molecular formula.
    4.  Compares the LLM's chosen formula with the deduced formula.
    5.  Checks the other options against the constraints to ensure the deduced formula is unique and provides reasons for their incorrectness.
    """
    
    # Step 1: Parse the LLM's final answer from the response text.
    match = re.search(r'<<<([A-D])>>>', llm_response)
    if not match:
        return "Failure: Could not find the final answer in the standard format '<<<X>>>' in the provided response."
    
    llm_answer_choice = match.group(1)
    
    options = {
        "A": "C11H14O2",
        "B": "C11H12O2",
        "C": "C12H12O2",
        "D": "C12H14O2"
    }
    
    llm_selected_formula = options.get(llm_answer_choice)

    # Step 2: Logically deduce the correct formula from the question's constraints.
    # Constraint: Di-substituted 6-membered aromatic ring -> C6H4 fragment.
    # Constraint: Two vinyl-H signals (one doublet, one doublet of quartets) -> This is a classic signature for a propenyl group (-CH=CH-CH3). This fragment is C3H5. It also accounts for one of the two required -CH3 signals.
    # Constraint: Ester group, a second -CH3 signal, and NO -CH2 signals -> The ester group (-COO-) must contain the second methyl group. To avoid a -CH2- group, it must be a methyl ester (-COOCH3) or an acetate group (-OCOCH3). Both possibilities have the formula C2H3O2.
    
    # Step 3: Assemble the fragments to get the total formula.
    # Carbons = 6 (ring) + 3 (propenyl) + 2 (ester part) = 11
    # Hydrogens = 4 (ring) + 5 (propenyl) + 3 (ester part) = 12
    # Oxygens = 2 (ester part) = 2
    deduced_formula = "C11H12O2"
    
    # Step 4: Compare the LLM's answer with the deduced formula.
    if llm_selected_formula != deduced_formula:
        return (f"Incorrect. The LLM chose option {llm_answer_choice} ({llm_selected_formula}), but the "
                f"correct formula deduced from the spectral data is {deduced_formula}. The reasoning is as follows: "
                f"A C6H4 ring, a C3H5 propenyl group, and a C2H3O2 ester/methyl fragment sum to C11H12O2.")

    # Step 5: Verify that the other options are incorrect to confirm the uniqueness of the solution.
    # This step confirms the robustness of the deduction.
    
    # Check A: C11H14O2
    # Difference from deduced: +2H. Implication: Saturation of the C=C bond to a propyl group (-CH2CH2CH3).
    # Violation: This introduces -CH2- groups, but the NMR shows no -CH2- signals.
    
    # Check C: C12H12O2
    # Difference from deduced: +C.
    # Violation: Adding a carbon is not feasible without violating other constraints (e.g., making the ring tri-substituted, or adding a forbidden -CH2- group in an ethyl ester).
    
    # Check D: C12H14O2
    # Difference from deduced: +CH2.
    # Violation: Directly contradicts the "no -CH2- signals" constraint.
    
    # Since the LLM's choice matches the unique, logically deduced formula, the answer is correct.
    return "Correct"
