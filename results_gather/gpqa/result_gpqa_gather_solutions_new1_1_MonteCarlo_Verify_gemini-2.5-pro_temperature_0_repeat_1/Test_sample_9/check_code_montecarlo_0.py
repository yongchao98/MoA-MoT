import math

def check_exoplanet_density_answer():
    """
    Checks the correctness of the answer about exoplanet density.

    The core principle is that for planets of the same composition, density
    increases with mass due to gravitational compression. This can be modeled
    with a scaling law.
    """

    # --- Define Constants and Model Parameters ---

    # Earth's average density in g/cm^3, our baseline.
    rho_earth = 5.51

    # For rocky planets, the radius (R) scales with mass (M) as R ∝ M^k.
    # The exponent 'k' is empirically and theoretically found to be around 0.27-0.3.
    # If k=1/3, density would be constant. Since k < 1/3, density increases with mass.
    # We'll use a standard value for the model.
    mass_radius_exponent_k = 0.27

    # From R ∝ M^k and ρ = M/V ∝ M/R^3, we get ρ ∝ M / (M^k)^3 = M^(1-3k).
    density_scaling_exponent = 1 - 3 * mass_radius_exponent_k

    # --- Calculate/Define Densities for Each Planet Option ---
    densities = {}

    # a) An Earth-mass and Earth-radius planet.
    # By definition, its density is Earth's density.
    densities['a'] = rho_earth

    # b) A planet with 2 Earth masses and a density of approximately 5.5 g/cm^3.
    # The density is explicitly given.
    densities['b'] = 5.5

    # c) A planet with the same composition as Earth but 5 times more massive.
    # We apply the scaling law: ρ_new = ρ_earth * (M_new / M_earth)^exponent
    mass_c_earth_masses = 5.0
    densities['c'] = rho_earth * (mass_c_earth_masses ** density_scaling_exponent)

    # d) A planet with the same composition as Earth but half the mass of Earth.
    mass_d_earth_masses = 0.5
    densities['d'] = rho_earth * (mass_d_earth_masses ** density_scaling_exponent)

    # --- Verify the Answer ---

    # Find which planet description ('a', 'b', 'c', or 'd') corresponds to the highest density.
    calculated_densest_planet = max(densities, key=densities.get)

    # The provided final answer is 'A'.
    provided_answer_choice = 'A'
    
    # Map the final answer choice to the planet description.
    # A) c, B) b, C) a, D) d
    answer_mapping = {'A': 'c', 'B': 'b', 'C': 'a', 'D': 'd'}
    claimed_densest_planet = answer_mapping.get(provided_answer_choice)

    # Check if the calculated result matches the claim from the provided answer.
    if calculated_densest_planet == claimed_densest_planet:
        return "Correct"
    else:
        reason = (
            f"Incorrect. The provided answer is '{provided_answer_choice}', which corresponds to planet '{claimed_densest_planet}'.\n"
            f"However, the physical model shows that planet '{calculated_densest_planet}' has the highest density.\n"
            f"This is because a more massive planet with the same composition experiences greater gravitational compression.\n\n"
            f"Calculated Densities (g/cm^3):\n"
            f"  - Planet a: {densities['a']:.2f}\n"
            f"  - Planet b: {densities['b']:.2f}\n"
            f"  - Planet c: {densities['c']:.2f} (Highest)\n"
            f"  - Planet d: {densities['d']:.2f}"
        )
        return reason

# Run the check and print the result.
result = check_exoplanet_density_answer()
print(result)