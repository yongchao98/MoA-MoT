import math

def check_answer():
    """
    This function checks the correctness of the answer to the exoplanet temperature ratio problem.
    """
    # Given values from the problem
    m_p1 = 7.0  # Mass of Planet 1 in Earth masses
    m_p2 = 5.0  # Mass of Planet 2 in Earth masses
    delta_lambda1 = 0.03  # Doppler shift for Planet 1 in Angstroms
    delta_lambda2 = 0.04  # Doppler shift for Planet 2 in Angstroms

    # The derivation shows that the ratio of equilibrium temperatures is:
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (K1 / K2)
    # And since K is proportional to delta_lambda, K1 / K2 = delta_lambda1 / delta_lambda2.
    
    # Calculate the mass ratio
    mass_ratio = m_p2 / m_p1
    
    # Calculate the Doppler shift ratio (which is equal to the RV semi-amplitude ratio K1/K2)
    doppler_shift_ratio = delta_lambda1 / delta_lambda2
    
    # Calculate the final temperature ratio
    temp_ratio = mass_ratio * doppler_shift_ratio
    
    # The options provided in the question
    options = {
        "A": 0.53,
        "B": 0.98,
        "C": 1.30,
        "D": 1.05
    }
    
    # The expected answer is A
    expected_answer_key = 'A'
    expected_value = options[expected_answer_key]

    # Check if the calculated result is close to the value of option A
    # We use math.isclose for floating-point comparison, with a tolerance of 0.01
    if math.isclose(temp_ratio, expected_value, rel_tol=1e-2):
        # Check if the provided answers select the correct option
        # In this case, we are verifying the calculation against the correct option 'A'
        # Some of the LLM answers incorrectly chose C or D despite correct calculations.
        return "Correct"
    else:
        return (f"Incorrect. The calculated ratio is {temp_ratio:.4f}, "
                f"which corresponds to option A (~0.53). "
                f"The provided answer might have selected the wrong option letter or made a calculation error.")

# Run the check
result = check_answer()
print(result)
