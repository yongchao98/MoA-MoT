def check_answer_correctness():
    """
    Checks the correctness of the answer to the textile dye question.

    The function verifies the answer based on the physical principle of fluorescence
    and the Stokes shift, which dictates that absorbed light must have a higher
    energy (shorter wavelength) than the emitted light.
    """
    # --- Problem & Answer Data ---
    energy_emitted_ev = 2.3393
    options = {
        'A': 'Violet',
        'B': 'Red',
        'C': 'Yellow',
        'D': 'Blue'
    }
    # The final answer provided by the LLM to be checked.
    llm_answer_letter = 'D'

    # --- Scientific Constants & Data ---
    # Approximate wavelength ranges for visible light in nanometers (nm).
    color_wavelengths = {
        'Violet': (380, 450),
        'Blue': (450, 495),
        'Green': (495, 570),
        'Yellow': (570, 590),
        'Red': (620, 750)
    }
    # Order of colors by energy (high to low)
    energy_order = ['Violet', 'Blue', 'Green', 'Yellow', 'Red']

    # --- Step 1: Calculate the wavelength of the emitted light ---
    # Using the formula: λ (nm) ≈ 1240 / E (eV)
    wavelength_emitted_nm = 1240 / energy_emitted_ev

    # --- Step 2: Determine the color of the emitted light ---
    emitted_color = None
    for color, (low, high) in color_wavelengths.items():
        if low <= wavelength_emitted_nm < high:
            emitted_color = color
            break

    if emitted_color != 'Green':
        return (f"Reasoning Error: The calculated emitted wavelength is {wavelength_emitted_nm:.2f} nm, "
                f"which corresponds to {emitted_color}, not Green. This affects the subsequent logic.")

    # --- Step 3: Apply the Stokes Shift principle ---
    # The absorbed light must have a shorter wavelength than the emitted light.
    # λ_absorbed < 530 nm
    possible_absorbed_colors = []
    for color, (low, high) in color_wavelengths.items():
        # The entire range of the color's wavelength must be shorter than the emitted wavelength.
        if high < wavelength_emitted_nm:
            possible_absorbed_colors.append(color)

    # --- Step 4: Determine the most plausible answer ---
    # In fluorescence, absorption and emission spectra are often adjacent.
    # The most likely absorbed color is the one with the next highest energy to the emitted color.
    try:
        emitted_color_index = energy_order.index(emitted_color)
        # The color before 'Green' in the energy_order list is 'Blue'.
        most_plausible_color = energy_order[emitted_color_index - 1]
    except (ValueError, IndexError):
        return "Internal Error: Could not determine the most plausible color from the energy order."

    # --- Step 5: Check the LLM's answer ---
    llm_chosen_color = options.get(llm_answer_letter)

    if llm_chosen_color == most_plausible_color:
        return "Correct"
    else:
        # Provide a specific reason for the incorrectness.
        if llm_chosen_color in ['Red', 'Yellow']:
            return (f"Incorrect. The answer '{llm_chosen_color}' violates the Stokes shift principle. "
                    f"The dye emits green light (~{wavelength_emitted_nm:.2f} nm), so it must absorb light with a "
                    f"SHORTER wavelength (higher energy). '{llm_chosen_color}' has a LONGER wavelength.")
        elif llm_chosen_color == 'Violet':
            return (f"Incorrect. While '{llm_chosen_color}' is a physically possible answer (shorter wavelength), "
                    f"it is not the most plausible one. For fluorescent dyes, the absorption spectrum is typically "
                    f"adjacent to the emission spectrum. The color with the next highest energy to the emitted 'Green' light is 'Blue'.")
        elif llm_chosen_color == 'Red' and emitted_color == 'Green':
             return (f"Incorrect. The answer '{llm_chosen_color}' is derived from the wrong physical principle. "
                     f"It incorrectly uses the 'complementary color' theory for reflected light. For emitted light (fluorescence), "
                     f"the Stokes shift principle must be used.")
        else:
            return f"Incorrect. The provided answer '{llm_chosen_color}' is not the correct choice."

# Run the check
result = check_answer_correctness()
print(result)