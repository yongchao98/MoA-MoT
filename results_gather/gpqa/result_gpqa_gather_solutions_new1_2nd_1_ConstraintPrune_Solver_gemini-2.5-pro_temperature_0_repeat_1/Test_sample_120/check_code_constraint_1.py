def check_chemistry_answer():
    """
    This function checks the correctness of the answer to the given organic chemistry question.
    It encodes the rules of the reaction (regioselectivity and stereoselectivity) to
    determine the expected product and compares it against the provided candidate answer.
    """

    # --- Problem Definition ---
    # Candidate answers provided in the question.
    options = {
        'A': {
            'name': '(1R,4R,5R)-2,2,4,5-tetramethylcyclohexan-1-ol',
            'structure': '2,2,4,5-tetramethylcyclohexan-1-ol',
            'stereochem': {'1': 'R', '4': 'R', '5': 'R'}
        },
        'B': {
            'name': '(1S,4R,5S)-2,2,4,5-tetramethylcyclohexan-1-ol',
            'structure': '2,2,4,5-tetramethylcyclohexan-1-ol',
            'stereochem': {'1': 'S', '4': 'R', '5': 'S'}
        },
        'C': {
            'name': '(1R,2S,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol',
            'structure': '1,2,4,5-tetramethylcyclohexan-1-ol',
            'stereochem': {'1': 'R', '2': 'S', '4': 'R', '5': 'R'}
        },
        'D': {
            'name': '(1R,2R,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol',
            'structure': '1,2,4,5-tetramethylcyclohexan-1-ol',
            'stereochem': {'1': 'R', '2': 'R', '4': 'R', '5': 'R'}
        }
    }

    # The final answer from the LLM to be checked.
    answer_to_check = 'C'

    # --- Constraint 1: Regioselectivity Analysis ---
    # The rule is that the nucleophile (Me-) attacks the less hindered carbon of the epoxide.
    # C1 is quaternary (more hindered), while C6 is tertiary (less hindered).
    # Attack occurs at C6, leading to a '1,2,4,5-tetramethylcyclohexan-1-ol' skeleton.
    expected_structure = '1,2,4,5-tetramethylcyclohexan-1-ol'

    candidate_structure = options[answer_to_check]['structure']
    if candidate_structure != expected_structure:
        return (f"Incorrect: The answer '{answer_to_check}' violates the regioselectivity constraint. "
                f"The nucleophile attacks the less hindered carbon (C6), which leads to a '{expected_structure}' skeleton, "
                f"but the answer proposes a '{candidate_structure}' skeleton, which would result from an incorrect attack at C1.")

    # --- Constraint 2: Stereoselectivity Analysis ---
    # Rule 1: Unaffected stereocenters (C1, C3, C4) retain their configuration.
    # Original C1(R) -> New C1(R)
    # Original C3(R) -> New C5(R)
    # Original C4(R) -> New C4(R)
    # Rule 2: The attacked stereocenter (C6) undergoes inversion.
    # Original C6 has (S) configuration. Due to geometric inversion and a change in CIP priorities
    # (high-priority epoxide oxygen replaced by lower-priority methyl), the final R/S descriptor is (S).
    expected_stereochem = {'1': 'R', '2': 'S', '4': 'R', '5': 'R'}

    candidate_stereochem = options[answer_to_check]['stereochem']

    if candidate_stereochem != expected_stereochem:
        for center, config in expected_stereochem.items():
            if candidate_stereochem.get(center) != config:
                return (f"Incorrect: The answer '{answer_to_check}' violates the stereoselectivity constraint. "
                        f"The configuration at C{center} should be ({config}) due to reaction rules, "
                        f"but the answer proposes ({candidate_stereochem.get(center)}).")
        # Fallback error message
        return (f"Incorrect: The answer '{answer_to_check}' has the wrong stereochemistry. "
                f"Expected {expected_stereochem} but got {candidate_stereochem}.")

    # --- Conclusion ---
    # If both constraints are satisfied, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_chemistry_answer()
print(result)