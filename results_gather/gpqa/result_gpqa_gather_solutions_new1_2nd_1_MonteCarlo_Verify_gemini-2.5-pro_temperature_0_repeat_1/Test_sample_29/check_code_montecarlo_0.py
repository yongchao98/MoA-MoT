import math

def check_answer():
    """
    Checks the correctness of the final answer by verifying each step of the physical reasoning.
    """
    # Step 0: Define the problem's parameters and the given answer.
    # The question provides a fixed mapping for the options.
    options = {
        "A": {"color": "Red", "min_wl_nm": 620, "max_wl_nm": 750},
        "B": {"color": "Yellow", "min_wl_nm": 570, "max_wl_nm": 590},
        "C": {"color": "Violet", "min_wl_nm": 380, "max_wl_nm": 450},
        "D": {"color": "Blue", "min_wl_nm": 450, "max_wl_nm": 495},
    }
    energy_emitted_eV = 2.3393
    final_answer_letter = "D"

    # Step 1: Calculate the wavelength of the emitted light.
    # The formula is λ (nm) ≈ 1240 / E (eV).
    wavelength_emitted_nm = 1240 / energy_emitted_eV
    
    # This wavelength (~530 nm) corresponds to green light.
    # We can add a check for this calculation.
    if not math.isclose(wavelength_emitted_nm, 530.07, rel_tol=1e-2):
        return f"Incorrect Calculation: The first step of calculating the emitted wavelength is wrong. Expected ~530.07 nm, but got {wavelength_emitted_nm:.2f} nm."

    # Step 2: Apply the Stokes Shift principle.
    # For fluorescence, the absorbed wavelength must be SHORTER than the emitted wavelength.
    # λ_absorbed < wavelength_emitted_nm
    
    valid_candidates = {}
    for letter, data in options.items():
        # The entire wavelength range of a valid candidate must be shorter than the emitted wavelength.
        if data["max_wl_nm"] < wavelength_emitted_nm:
            valid_candidates[letter] = data

    # Step 3: Check if the provided answer is among the valid candidates.
    if final_answer_letter not in valid_candidates:
        chosen_data = options[final_answer_letter]
        return (f"Incorrect: The final answer '{chosen_data['color']}' ({final_answer_letter}) violates the Stokes Shift principle. "
                f"Its wavelength range ({chosen_data['min_wl_nm']}-{chosen_data['max_wl_nm']} nm) is not shorter than the emitted wavelength of {wavelength_emitted_nm:.2f} nm (Green). "
                "The absorbed light must have higher energy (shorter wavelength) than the emitted light.")

    # Step 4: Apply the plausibility rule to find the *best* answer.
    # The emitted color is Green. The valid absorbed colors are Blue and Violet.
    # The most plausible absorbed color is the one adjacent to the emitted color on the higher-energy side of the spectrum.
    # The order of colors by increasing energy is: Red < Yellow < Green < Blue < Violet.
    # The color adjacent to Green with higher energy is Blue.
    most_plausible_color = "Blue"
    
    chosen_color = options[final_answer_letter]["color"]
    
    if chosen_color != most_plausible_color:
        return (f"Incorrect: While '{chosen_color}' is an energetically possible answer, it is not the most plausible one. "
                f"The emitted light is green (~530 nm). The most plausible absorbed color is the adjacent one with higher energy, which is '{most_plausible_color}', not '{chosen_color}'.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_answer()
print(result)