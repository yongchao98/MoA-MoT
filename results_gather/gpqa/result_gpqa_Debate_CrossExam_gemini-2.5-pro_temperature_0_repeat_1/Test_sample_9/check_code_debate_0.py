import math

def check_exoplanet_density_answer():
    """
    This function checks the correctness of the provided answer about exoplanet density.

    The question asks to identify the exoplanet with the highest density among four options:
    a) An Earth-mass and Earth-radius planet.
    b) A planet with 2 Earth masses and a density of approximately 5.5 g/cm^3.
    c) A planet with the same composition as Earth but 5 times more massive than Earth.
    d) A planet with the same composition as Earth but half the mass of Earth.

    The provided answer is 'D', which corresponds to option 'c'. This implies that
    planet 'c' should have the highest density.
    """

    # --- Model and Assumptions ---
    # The core principle is that for planets of the same composition, a more massive
    # planet will be denser due to stronger gravitational compression. We can quantify
    # this using a standard mass-radius relationship for rocky planets.
    #
    # A common approximation is R/R_earth = (M/M_earth)^β. For rocky planets, β is
    # less than 1/3 (the value for constant density). A typical value from models
    # (e.g., Zeng et al., 2016) is β ≈ 0.274.
    #
    # Density (ρ) is Mass (M) / Volume (V), and V ∝ R³.
    # The relative density ρ/ρ_earth can be expressed as:
    # ρ_rel = (M/M_earth) / (R/R_earth)³
    # ρ_rel = (M_rel) / (M_rel^β)³ = M_rel / M_rel^(3β) = M_rel^(1 - 3β)
    #
    # Using β = 0.274, the exponent for the density relationship is:
    # exponent = 1 - 3 * 0.274 = 1 - 0.822 = 0.178
    # So, our model is: ρ = ρ_earth * (M/M_earth)^0.178

    # --- Constants ---
    # Average density of Earth in g/cm^3
    rho_earth = 5.51

    # --- Calculate estimated densities for each option ---
    densities = {}

    # Option a: An Earth-mass and Earth-radius planet. This is our baseline.
    # Mass relative to Earth = 1.0
    densities['a'] = rho_earth * (1.0 ** 0.178)

    # Option b: A planet with 2 Earth masses and a density of approximately 5.5 g/cm^3.
    # The density is given directly in the problem statement.
    densities['b'] = 5.50

    # Option c: A planet with the same composition as Earth but 5 times more massive.
    # Mass relative to Earth = 5.0
    densities['c'] = rho_earth * (5.0 ** 0.178)

    # Option d: A planet with the same composition as Earth but half the mass.
    # Mass relative to Earth = 0.5
    densities['d'] = rho_earth * (0.5 ** 0.178)

    # --- Analysis ---
    # Find the option with the highest calculated density.
    try:
        highest_density_option = max(densities, key=densities.get)
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # The provided answer is 'D', which corresponds to option 'c'.
    # We check if our calculated densest planet matches the one from the answer.
    expected_densest_option = 'c'

    if highest_density_option == expected_densest_option:
        # The calculation confirms the reasoning. Let's double-check that all
        # constraints are met, i.e., that 'c' is denser than all other options.
        if densities['c'] > densities['a'] and densities['c'] > densities['b'] and densities['c'] > densities['d']:
             return "Correct"
        else:
            # This case is unlikely but represents a logical check.
            reason = "The qualitative reasoning is correct, but the quantitative check fails.\n"
            reason += f"Calculated densities (g/cm^3): a={densities['a']:.2f}, b={densities['b']:.2f}, c={densities['c']:.2f}, d={densities['d']:.2f}.\n"
            reason += f"Although option 'c' was expected to be the densest, the comparison shows this is not the case."
            return reason
    else:
        reason = f"The provided answer is incorrect.\n"
        reason += f"The analysis is based on the principle that for a fixed composition, density increases with mass due to gravitational compression.\n"
        reason += f"Calculated densities (g/cm^3):\n"
        reason += f"  - a) Earth-mass planet: {densities['a']:.2f}\n"
        reason += f"  - b) 2 Earth-mass planet: {densities['b']:.2f} (given)\n"
        reason += f"  - c) 5 Earth-mass planet: {densities['c']:.2f}\n"
        reason += f"  - d) 0.5 Earth-mass planet: {densities['d']:.2f}\n"
        reason += f"The planet with the highest calculated density is option '{highest_density_option}' with a density of {densities[highest_density_option]:.2f} g/cm^3.\n"
        reason += f"The provided answer states that option '{expected_densest_option}' is the densest, which contradicts the calculation."
        return reason

# The final result of the check.
result = check_exoplanet_density_answer()
print(result)