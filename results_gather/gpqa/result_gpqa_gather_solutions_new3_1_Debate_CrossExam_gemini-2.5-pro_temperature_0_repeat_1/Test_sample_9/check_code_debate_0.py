import math

def check_exoplanet_density():
    """
    Checks the correctness of the answer by modeling planetary densities.

    The model uses the principle of gravitational compression, where for planets
    of the same composition, density increases with mass. A common empirical
    mass-radius relationship for rocky planets is R ∝ M^β, where β is an
    exponent less than 1/3 (the value for a constant-density object).
    From this, density ρ ∝ M / R^3 ∝ M / (M^β)^3 ∝ M^(1 - 3β).
    Since β < 1/3, the exponent (1 - 3β) is positive, confirming that
    density increases with mass.
    """

    # --- 1. Define Model Parameters and Constants ---
    # Earth's average density in g/cm^3 as a baseline.
    rho_earth = 5.51

    # A typical mass-radius exponent (beta) for rocky planets. Any value
    # between 0 and 1/3 will yield the same qualitative result.
    # We use a common approximation.
    beta = 0.27

    # --- 2. Calculate the Density for Each Planet Option ---
    # a) An Earth-mass and Earth-radius planet.
    rho_a = rho_earth

    # b) A planet with 2 Earth masses and a density of approximately 5.5 g/cm^3.
    rho_b = 5.5

    # c) A planet with the same composition as Earth but 5 times more massive.
    # Its density relative to Earth is scaled by (Mass_Ratio)^(1 - 3*beta).
    mass_ratio_c = 5.0
    rho_c = rho_earth * math.pow(mass_ratio_c, 1 - 3 * beta)

    # d) A planet with the same composition as Earth but half the mass of Earth.
    mass_ratio_d = 0.5
    rho_d = rho_earth * math.pow(mass_ratio_d, 1 - 3 * beta)

    densities = {
        'a': rho_a,
        'b': rho_b,
        'c': rho_c,
        'd': rho_d
    }

    # --- 3. Identify the Densest Planet ---
    # Find the key ('a', 'b', 'c', or 'd') corresponding to the max value.
    densest_planet_description = max(densities, key=densities.get)

    # --- 4. Map the Description to the Final Answer Choice ---
    # The question maps descriptions to final choices as: A)a, B)d, C)b, D)c
    answer_mapping = {
        'a': 'A',
        'b': 'C',
        'c': 'D',
        'd': 'B'
    }
    calculated_correct_choice = answer_mapping[densest_planet_description]

    # The final answer provided by the LLM.
    llm_final_answer = "D"

    # --- 5. Verify the LLM's Answer ---
    if calculated_correct_choice == llm_final_answer:
        # The reasoning and the final answer are both correct.
        return "Correct"
    else:
        # The reasoning or the final mapping was incorrect.
        error_reason = (
            f"The final answer is incorrect. "
            f"The analysis shows that planet (c) is the densest. "
            f"Calculated densities (g/cm^3): a={rho_a:.2f}, b={rho_b:.2f}, c={rho_c:.2f}, d={rho_d:.2f}. "
            f"The densest planet is '{densest_planet_description}', which corresponds to answer choice '{calculated_correct_choice}'. "
            f"The provided answer was '{llm_final_answer}'."
        )
        return error_reason

# Execute the check and print the result.
result = check_exoplanet_density()
print(result)