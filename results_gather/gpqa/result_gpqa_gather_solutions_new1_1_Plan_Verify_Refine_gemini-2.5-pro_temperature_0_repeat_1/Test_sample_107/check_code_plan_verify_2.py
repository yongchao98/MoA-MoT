import scipy.constants as const
import math

def check_correctness_of_llm_answer():
    """
    This function checks the correctness of the LLM's answer to the physics problem.

    The problem asks to compare the order of magnitude of:
    1. Paramagnetic coupling term (<H>) for B=1T, small m.
    2. Transition energy (ΔE) for λ=0.4861μm.

    The possible answers are:
    A) <H> > ΔE
    B) <H> = ΔE
    C) <H> << ΔE
    D) <H> >> ΔE

    The LLM's final answer is C. This code will verify if C is the correct choice.
    """

    # --- Define problem parameters from the question ---
    B = 1.0  # Magnetic field in Tesla
    m = 1    # A representative small integer for the orbital magnetic quantum number
    lambda_val = 0.4861e-6  # Wavelength in meters

    # --- Use high-precision physical constants from scipy ---
    # Bohr magneton in J/T
    mu_B = const.physical_constants['Bohr magneton'][0]
    # Planck's constant in J*s
    h = const.h
    # Speed of light in m/s
    c = const.c

    # --- Perform the calculations ---
    # 1. Calculate the paramagnetic coupling energy <H>
    # Formula: <H> = m * μ_B * B
    H_coupling = m * mu_B * B

    # 2. Calculate the transition energy ΔE
    # Formula: ΔE = h * c / λ
    delta_E = (h * c) / lambda_val

    # --- Compare the two energies to determine the correct relationship ---
    # The most robust way to compare orders of magnitude is to calculate the ratio.
    if delta_E == 0:
        # Avoid division by zero, though physically impossible here.
        return "Error: Transition energy calculated to be zero."

    ratio = H_coupling / delta_E

    # Determine the correct option based on the ratio.
    # We define "much less" (<<) or "much greater" (>>) as a difference of at least
    # 2-3 orders of magnitude (a factor of 100-1000).
    
    determined_option = None
    if ratio < 1e-3:
        determined_option = "C"  # <H> << ΔE
    elif ratio > 1e3:
        determined_option = "D"  # <H> >> ΔE
    elif math.isclose(ratio, 1.0, rel_tol=0.1): # A generous tolerance for equality
        determined_option = "B"  # <H> = ΔE
    elif ratio < 1.0:
        # This case would be '<', but since '<<' (C) is an option and the ratio is
        # extremely small, it falls under C. The only remaining option is A.
        determined_option = "A"
    else: # ratio > 1.0
        determined_option = "A"  # <H> > ΔE

    # --- Verify against the LLM's answer ---
    # The provided answer from the meta-analysis is 'C'.
    llm_answer = "C"

    if determined_option == llm_answer:
        return "Correct"
    else:
        reason = (
            f"The answer is incorrect.\n"
            f"The calculation shows the correct option is {determined_option}, but the provided answer is {llm_answer}.\n"
            f"--- Calculation Details ---\n"
            f"Paramagnetic Coupling Energy (<H>): {H_coupling:.4e} J\n"
            f"Transition Energy (ΔE): {delta_E:.4e} J\n"
            f"Ratio (<H> / ΔE): {ratio:.4e}\n"
            f"Since the ratio is extremely small ({ratio:.2e}), the relationship is <H> << ΔE, which corresponds to option {determined_option}."
        )
        return reason

# Run the check and print the result
result = check_correctness_of_llm_answer()
print(result)