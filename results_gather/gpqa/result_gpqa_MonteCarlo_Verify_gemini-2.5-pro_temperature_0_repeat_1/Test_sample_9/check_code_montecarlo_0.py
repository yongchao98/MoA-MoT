import math

def check_exoplanet_density_answer():
    """
    Checks the correctness of the answer to the exoplanet density question.

    The question asks which of the following has the highest density:
    a) An Earth-mass and Earth-radius planet.
    b) A planet with 2 Earth masses and a density of approximately 5.5 g/cm^3.
    c) A planet with the same composition as Earth but 5 times more massive than Earth.
    d) A planet with the same composition as Earth but half the mass of Earth.

    The provided answer from the other LLM is 'C'. This checker verifies if option 'c'
    is indeed the one with the highest density.
    """

    # --- Define physical parameters and relationships ---

    # Earth's density in g/cm^3. This is our baseline for options a, c, and d.
    DENSITY_EARTH = 5.51

    # For planets of the same composition, density (rho) scales with mass (M)
    # according to the formula: rho_new = rho_ref * (M_new / M_ref)^(1 - 3*beta).
    # 'beta' is the mass-radius exponent (R ∝ M^beta). For rocky planets, beta is
    # empirically found to be less than 1/3 (typically ~0.27).
    # This means the exponent (1 - 3*beta) is a positive number, causing density
    # to increase with mass. We use a representative value for beta.
    BETA = 0.27
    exponent = 1 - 3 * BETA

    # --- Calculate the density for each option ---

    # Option a: Baseline Earth
    density_a = DENSITY_EARTH

    # Option b: Given density
    # The question states "approximately 5.5 g/cm^3".
    density_b = 5.5

    # Option c: 5x Earth mass, same composition
    mass_ratio_c = 5.0
    density_c = DENSITY_EARTH * (mass_ratio_c ** exponent)

    # Option d: 0.5x Earth mass, same composition
    mass_ratio_d = 0.5
    density_d = DENSITY_EARTH * (mass_ratio_d ** exponent)

    # --- Verify the answer ---

    # Store densities in a dictionary for easy comparison
    densities = {
        "a": density_a,
        "b": density_b,
        "c": density_c,
        "d": density_d,
    }

    # Find the option with the maximum calculated density
    calculated_highest_option = max(densities, key=densities.get)

    # The provided answer is 'C', which corresponds to option 'c'.
    provided_answer_option = 'c'

    if calculated_highest_option == provided_answer_option:
        # The calculation confirms the provided answer.
        # We can perform a final sanity check on the ordering.
        # The expected order is d < b < a < c.
        # Our calculated values are:
        # d ≈ 4.8 g/cm^3
        # b = 5.5 g/cm^3
        # a = 5.51 g/cm^3
        # c ≈ 7.5 g/cm^3
        # The order holds.
        return "Correct"
    else:
        # This block would execute if the calculation contradicted the provided answer.
        reason = (f"The answer is incorrect. "
                  f"Based on the physical model of gravitational compression, the planet with the highest density is option '{calculated_highest_option}'. "
                  f"The calculated densities are: a) {density_a:.2f}, b) {density_b:.2f}, c) {density_c:.2f}, d) {density_d:.2f} g/cm^3. "
                  f"The provided answer was 'C', but option '{calculated_highest_option}' has the highest density.")
        return reason

# Print the result of the check
print(check_exoplanet_density_answer())