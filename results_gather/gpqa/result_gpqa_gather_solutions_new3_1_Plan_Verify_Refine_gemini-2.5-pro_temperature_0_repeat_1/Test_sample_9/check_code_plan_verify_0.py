import math

def check_correctness():
    """
    Checks the correctness of the answer to the exoplanet density question.

    The core logic is based on the physical principle of gravitational compression.
    For planets of the same composition, density increases with mass. A common
    approximation for the mass-radius relationship is R ∝ M^β, where β ≈ 0.27.
    From this, we can derive the density relationship: ρ ∝ M^(1-3β).
    """
    # The final answer provided in the prompt to be checked.
    provided_answer = "C"

    # --- Step 1: Model the density of each planet ---

    # Baseline Earth density in g/cm^3
    earth_density = 5.51

    # Exponent for the mass-radius relationship (R ∝ M^beta) for rocky planets.
    # A value of 0.27 is a standard approximation.
    beta = 0.27

    densities = {}

    # Option a) An Earth-mass and Earth-radius planet. This is our baseline.
    densities['a'] = earth_density

    # Option b) A planet with 2 Earth masses and a density of approximately 5.5 g/cm^3.
    densities['b'] = 5.5

    # Option c) A planet with the same composition as Earth but 5 times more massive.
    # We calculate its density relative to Earth's using the formula: ρ_new = ρ_earth * (M_new/M_earth)^(1 - 3*β)
    mass_ratio_c = 5.0
    density_c = earth_density * (mass_ratio_c ** (1 - 3 * beta))
    densities['c'] = density_c

    # Option d) A planet with the same composition as Earth but half the mass of Earth.
    mass_ratio_d = 0.5
    density_d = earth_density * (mass_ratio_d ** (1 - 3 * beta))
    densities['d'] = density_d

    # --- Step 2: Determine the correct answer based on the model ---

    # Find the description ('a', 'b', 'c', or 'd') corresponding to the highest density.
    if not densities:
        return "Error: Densities dictionary is empty."
        
    highest_density_description = max(densities, key=densities.get)

    # Map the correct description to the final answer choice (A, B, C, D).
    # A) b, B) d, C) c, D) a
    description_to_choice_map = {
        'a': 'D',
        'b': 'A',
        'c': 'C',
        'd': 'B'
    }
    
    correct_choice = description_to_choice_map.get(highest_density_description)

    # --- Step 3: Compare with the provided answer and return the result ---

    if correct_choice == provided_answer:
        return "Correct"
    else:
        reason = (
            f"The provided answer '{provided_answer}' is incorrect.\n"
            f"The analysis shows that planet '{highest_density_description}' has the highest density.\n"
            f"Calculated densities:\n"
            f"  - a: {densities['a']:.2f} g/cm^3\n"
            f"  - b: {densities['b']:.2f} g/cm^3\n"
            f"  - c: {densities['c']:.2f} g/cm^3\n"
            f"  - d: {densities['d']:.2f} g/cm^3\n"
            f"The highest density belongs to planet '{highest_density_description}', which corresponds to choice '{correct_choice}'. "
            f"Therefore, the correct answer is '{correct_choice}'."
        )
        return reason

# Execute the check and print the result.
print(check_correctness())