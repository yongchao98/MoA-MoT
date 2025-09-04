import math

def check_exoplanet_density():
    """
    This function evaluates the density of four hypothetical exoplanets
    to determine which one is the densest, thereby checking the problem's answer.
    """
    # --- Define Baseline and Given Values ---

    # Option a) An Earth-mass and Earth-radius planet.
    # This is Earth. Its average density is a known constant.
    # Source: NASA Earth Fact Sheet
    density_earth = 5.514  # in g/cm^3
    density_a = density_earth

    # Option b) A planet with 2 Earth masses and a density of approximately 5.5 g/cm^3.
    # The density is given directly in the problem.
    density_b = 5.5

    # --- Model for Options c and d ---

    # For planets of the same composition, density increases with mass due to
    # gravitational compression. We can model this with the mass-radius relationship
    # R ∝ M^k, where k is approx. 0.27 for rocky planets.
    # This means density ρ ∝ M^(1-3k), so density increases with mass.
    k = 0.27
    exponent = 1 - 3 * k

    # Option c) A planet with the same composition as Earth but 5 times more massive.
    mass_ratio_c = 5.0
    density_c = density_earth * (mass_ratio_c ** exponent)

    # Option d) A planet with the same composition as Earth but half the mass of Earth.
    mass_ratio_d = 0.5
    density_d = density_earth * (mass_ratio_d ** exponent)

    # --- Comparison and Conclusion ---

    densities = {
        'a': density_a,
        'b': density_b,
        'c': density_c,
        'd': density_d
    }

    print("--- Density Calculation ---")
    print(f"a) Earth-like planet: {densities['a']:.2f} g/cm^3")
    print(f"b) 2 M_earth planet: {densities['b']:.2f} g/cm^3 (given)")
    print(f"c) 5 M_earth planet (same composition): {densities['c']:.2f} g/cm^3 (calculated)")
    print(f"d) 0.5 M_earth planet (same composition): {densities['d']:.2f} g/cm^3 (calculated)")
    print("-" * 27)

    # Find the option with the highest density
    highest_density_option = max(densities, key=densities.get)
    highest_density_value = densities[highest_density_option]

    print(f"Conclusion: Option '{highest_density_option}' has the highest density at {highest_density_value:.2f} g/cm^3.")

    # This check confirms that 'c' is the correct answer.
    if highest_density_option == 'c':
        print("\nThe code confirms that 'C' is the correct answer.")
    else:
        print(f"\nThe code indicates '{highest_density_option.upper()}' is the correct answer, which contradicts the expected physical result.")

# Run the check
check_exoplanet_density()