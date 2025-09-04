import math

def check_potential_energy_answer():
    """
    This function checks the correctness of the provided LLM answer for the potential energy
    of a charge 'q' near a grounded conducting sphere of radius 'R'.

    The LLM's answer provides a derivation and concludes that option D is correct.
    Option D: U = -(1/2) * k*q^2*R / (d^2 - R^2)

    This code will:
    1.  Use a set of numerical values for the physical parameters (k, q, R, d).
    2.  Calculate the expected potential energy using the well-established formula from physics,
        which matches the derivation provided by the LLM.
    3.  Evaluate the formula for the chosen option (D).
    4.  Compare the two results to confirm they match.
    5.  Verify that other options (A, B, C) do not match, ensuring the answer is unambiguous.
    """
    # Let's use a set of simple, non-trivial values for the physical parameters.
    # We can set k=1 and q=1 without loss of generality to check the formula's structure.
    k = 1.0
    q = 1.0
    # Choose R and d such that the charge is outside the sphere (d > R).
    R = 2.0
    d = 4.0

    # The provided answer from the LLM is option D.
    llm_selected_option = 'D'

    # The correct formula for the potential energy of this system is derived
    # from the method of images. The energy is half the interaction energy
    # of the charge with its image.
    # U = (1/2) * q * V_image = -(1/2) * k * q^2 * R / (d^2 - R^2)
    try:
        ground_truth_value = -0.5 * k * q**2 * R / (d**2 - R**2)
    except ZeroDivisionError:
        return "Constraint violated: d cannot be equal to R."


    # Define the formulas for all given options.
    options = {
        'A': -k * q**2 * d / (d**2 - R**2),
        'B': -0.5 * k * q**2 * R**2 / (d**2 - R**2),
        'C': -0.5 * k * q**2 * d / (d**2 + R**2),
        'D': -0.5 * k * q**2 * R / (d**2 - R**2)
    }

    # Get the value calculated from the formula of the selected option.
    value_from_selected_option = options[llm_selected_option]

    # 1. Check if the selected option's value matches the ground truth.
    if not math.isclose(value_from_selected_option, ground_truth_value, rel_tol=1e-9):
        # This would mean the LLM's derivation or its choice of option is wrong.
        # Let's find which option (if any) is actually correct.
        correct_option_key = None
        for key, value in options.items():
            if math.isclose(value, ground_truth_value, rel_tol=1e-9):
                correct_option_key = key
                break
        
        if correct_option_key:
            return (f"Incorrect. The provided answer is '{llm_selected_option}', but the correct formula "
                    f"corresponds to option '{correct_option_key}'.")
        else:
            return ("Incorrect. The provided answer is '{llm_selected_option}', and none of the options "
                    "match the correct physical formula.")

    # 2. Verify that the correct formula is unique among the options.
    # This ensures the question is well-posed for the given parameters.
    for key, value in options.items():
        if key != llm_selected_option and math.isclose(value, ground_truth_value, rel_tol=1e-9):
            return (f"Ambiguous Answer. The selected option '{llm_selected_option}' is correct, but option '{key}' "
                    f"yields the same numerical result for the tested parameters. The formulas may be algebraically equivalent "
                    f"under certain conditions or the test case is degenerate.")

    # 3. If the selected option matches the ground truth and is unique, the answer is correct.
    return "Correct"

# Execute the check and print the result.
print(check_potential_energy_answer())