import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the provided answer to the exoplanet physics problem.

    The core physics principles are:
    1. The radial velocity semi-amplitude (K) is directly proportional to the observed wavelength shift (Δλ).
       Therefore, K₁ / K₂ = Δλ₁ / Δλ₂.
    2. For a circular orbit, K is proportional to P^(-1/3), where P is the orbital period, assuming
       star mass, planet mass, and orbital inclination are constant between the two systems.
    3. Combining these, we get (P₂ / P₁)^(1/3) = K₁ / K₂ = Δλ₁ / Δλ₂.
    4. Solving for the desired ratio gives P₂ / P₁ = (Δλ₁ / Δλ₂)^3.
    """

    # Given values from the question
    delta_lambda_1 = 5  # miliangstrom
    delta_lambda_2 = 7  # miliangstrom

    # The options as listed in the final answer's analysis
    options = {
        'A': 1.40,
        'B': 1.96,
        'C': 0.36,
        'D': 0.85
    }

    # The final answer provided by the agent
    agent_answer_choice = 'C'
    
    # --- Step 1: Calculate the correct numerical answer based on physics principles ---
    try:
        # The ratio of the periods P₂ / P₁ is (Δλ₁ / Δλ₂)^3
        correct_ratio = (delta_lambda_1 / delta_lambda_2) ** 3
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Step 2: Check if the agent's chosen option corresponds to the calculated value ---
    agent_answer_value = options.get(agent_answer_choice)

    if agent_answer_value is None:
        return f"The agent's final answer '{agent_answer_choice}' is not a valid option."

    # The options are approximate, so we check if the calculated value is close to the option's value.
    # A relative tolerance of 5% is reasonable for "approximately equal to".
    if not math.isclose(correct_ratio, agent_answer_value, rel_tol=0.05):
        # Find the option that is actually the closest to the correct answer
        closest_option = min(options.keys(), key=lambda k: abs(options[k] - correct_ratio))
        
        return (f"Incorrect. The calculated ratio P₂ / P₁ is (5/7)^3 ≈ {correct_ratio:.4f}. "
                f"This value is closest to option '{closest_option}' ({options[closest_option]}). "
                f"The agent chose option '{agent_answer_choice}' ({agent_answer_value}), which is not the best match.")

    # --- Step 3: Verify the reasoning provided by the agent ---
    # The agent's reasoning is: P₂ / P₁ = (5 / 7)^3 ≈ 0.3644, which matches option C (~0.36).
    # This reasoning is sound and matches our calculation.
    
    return "Correct"

# Run the check
result = check_correctness_of_answer()
print(result)