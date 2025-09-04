import sys
import io

def check_correctness():
    """
    This function checks the correctness of the LLM's answer by following these steps:
    1. Defines the physical constants and color spectrum data.
    2. Calculates the wavelength of the emitted light from the given energy.
    3. Verifies that the emitted light is green.
    4. Applies the Stokes Shift principle: the absorbed light must have a shorter wavelength (higher energy) than the emitted light.
    5. Filters the options to find all physically possible answers.
    6. Selects the most plausible answer from the valid options based on Kasha's rule (fluorescence typically originates from the lowest-energy absorption band).
    7. Compares this result with the provided answer.
    """
    # --- Step 1: Define constants and data ---
    # Physical constant for E(eV) to lambda(nm) conversion
    hc_approx = 1240  # eV*nm
    # Given energy of emitted light
    emitted_energy_eV = 2.3393
    
    # Approximate wavelength ranges for visible light (in nm), based on the consensus in the provided answers.
    color_wavelengths = {
        "Violet": (400, 450),
        "Blue": (450, 495),
        "Green": (495, 570),
        "Yellow": (570, 590),
        "Red": (620, 750)
    }
    
    # Question options
    options = {
        "A": "Violet",
        "B": "Yellow",
        "C": "Red",
        "D": "Blue"
    }
    
    # The final answer from the LLM to be checked
    llm_answer_key = "D"
    llm_answer_color = options[llm_answer_key]

    # --- Step 2: Calculate emitted wavelength ---
    emitted_wavelength_nm = hc_approx / emitted_energy_eV
    
    # --- Step 3: Verify emitted color ---
    emitted_color = None
    for color, (min_wl, max_wl) in color_wavelengths.items():
        if min_wl <= emitted_wavelength_nm <= max_wl:
            emitted_color = color
            break
            
    if emitted_color != "Green":
        return f"Reason: The initial calculation is flawed. The emitted light with wavelength {emitted_wavelength_nm:.2f} nm is {emitted_color}, not Green as reasoned."

    # --- Step 4: Apply Stokes Shift and find valid options ---
    # The absorbed wavelength must be shorter than the emitted wavelength.
    valid_options = {}
    for key, color in options.items():
        # The entire range of the absorbed color should be at a shorter wavelength.
        # We check if the maximum wavelength of the color is less than the emitted wavelength.
        _min_wl, max_wl = color_wavelengths[color]
        if max_wl < emitted_wavelength_nm:
            valid_options[key] = color
            
    # Check if the LLM's answer is among the physically possible options
    if llm_answer_key not in valid_options:
        return f"Reason: The answer '{llm_answer_color}' is incorrect because it violates the Stokes Shift principle. Its wavelength range ({color_wavelengths[llm_answer_color][0]}-{color_wavelengths[llm_answer_color][1]} nm) is not shorter than the emitted wavelength of {emitted_wavelength_nm:.2f} nm."

    # --- Step 5: Select the most plausible answer ---
    # Kasha's rule implies fluorescence occurs from the lowest excited state (S1).
    # Therefore, the most likely absorption is the lowest-energy (longest wavelength) one that is still higher in energy than the emission.
    # We find the valid option with the longest wavelength.
    most_plausible_color = None
    longest_wavelength_in_valid_options = 0
    
    for key, color in valid_options.items():
        _min_wl, max_wl = color_wavelengths[color]
        if max_wl > longest_wavelength_in_valid_options:
            longest_wavelength_in_valid_options = max_wl
            most_plausible_color = color
            
    # --- Step 6: Compare and return result ---
    if most_plausible_color == llm_answer_color:
        return "Correct"
    else:
        return f"Reason: The answer '{llm_answer_color}' is incorrect. While it is a physically possible option, '{most_plausible_color}' is the most plausible answer. For fluorescence, the absorption band is typically the one with the lowest energy (longest wavelength) that is still higher than the emission energy. Among the valid options {list(valid_options.values())}, '{most_plausible_color}' has the longest wavelength."

# Execute the check
result = check_correctness()
print(result)