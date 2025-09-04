import scipy.constants as const
import re

def check_physics_answer():
    """
    Checks the correctness of the LLM's answer to the physics problem.

    The function calculates the two energies in question:
    1. The transition energy (ΔE) for a photon of a given wavelength.
    2. The paramagnetic coupling energy (<H>) for an atom in a given magnetic field.

    It then compares their orders of magnitude to determine the correct relationship
    and verifies if the provided answer matches this conclusion.
    """
    # --- Given values from the question ---
    # Wavelength in meters
    wavelength = 0.4861e-6  # 0.4861 μm
    # Magnetic field strength in Tesla
    B_field = 1.0
    # Orbital magnetic quantum number (assumed to be a small integer, m=1)
    m_quantum_number = 1

    # --- Physical Constants from scipy.constants for high precision ---
    # Planck's constant in J·s
    h = const.h
    # Speed of light in m/s
    c = const.c
    # Bohr magneton in J/T
    mu_B = const.value('Bohr magneton')

    # --- Step 1: Calculate the transition energy (ΔE) ---
    # Formula: ΔE = hc/λ
    try:
        transition_energy = (h * c) / wavelength
    except ZeroDivisionError:
        return "Error: Wavelength cannot be zero."

    # --- Step 2: Calculate the paramagnetic coupling energy (<H>) ---
    # Formula: <H> = m * μ_B * B
    paramagnetic_energy = m_quantum_number * mu_B * B_field

    # --- Step 3: Compare the two energies ---
    # The most robust way to compare orders of magnitude is to find their ratio.
    if transition_energy == 0:
        return "Error: Calculated transition energy is zero, cannot compute ratio."
    ratio = paramagnetic_energy / transition_energy

    # --- Step 4: Determine the correct relationship based on the ratio ---
    # The options given in the prompt are:
    # A) <H> > ΔE
    # B) <H> = ΔE
    # C) <H> >> ΔE
    # D) <H> << ΔE
    # A ratio of ~10^-5 clearly indicates that <H> is "much, much less" than ΔE.
    # We set a threshold for "much less" (e.g., < 0.001) and "much greater" (e.g., > 1000).
    
    correct_option = None
    if ratio < 1e-3:
        correct_option = 'D'  # <H> << ΔE
    elif ratio > 1e3:
        correct_option = 'C'  # <H> >> ΔE
    elif 1 - 1e-3 < ratio < 1 + 1e-3:
        correct_option = 'B'  # <H> ≈ ΔE
    elif ratio > 1:
        correct_option = 'A'  # <H> > ΔE
    else:
        # This case would be for <H> < ΔE but not << ΔE, which is not an option.
        # Given the options, D is the only logical choice for ratio < 1.
        pass

    # --- Step 5: Check the provided answer ---
    # The final answer from the LLM's analysis is <<<D>>>
    llm_answer = "<<<D>>>"
    
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        return f"Invalid answer format provided: {llm_answer}"
    
    llm_choice = match.group(1)

    if llm_choice == correct_option:
        return "Correct"
    else:
        return (f"Incorrect. The calculation shows that the paramagnetic coupling energy "
                f"<H> is approximately {paramagnetic_energy:.3e} J, and the transition energy "
                f"ΔE is approximately {transition_energy:.3e} J. The ratio <H>/ΔE is {ratio:.3e}, "
                f"which means <H> << ΔE. This corresponds to option {correct_option}, but the provided answer was {llm_choice}.")

# Execute the check and print the result
result = check_physics_answer()
print(result)