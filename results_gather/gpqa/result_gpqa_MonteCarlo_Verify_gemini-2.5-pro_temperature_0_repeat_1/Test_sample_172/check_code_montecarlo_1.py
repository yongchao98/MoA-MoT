import math

def check_correctness():
    """
    This function checks the correctness of the given answer to the physics problem.
    
    Problem:
    If uncertainty in space of electron's location, which is travelling with speed v= 2* 10^8 m/s 
    along x-direction is Δx=0.1 nm . Based on the information estimate the minimum 
    uncertainty in the energy ΔE of electron.

    Options:
    A) ~10^(-17) J
    B) ~10^(-19) J
    C) ~10^(-16) J
    D) ~10^(-18) J

    Provided Answer from LLM: C
    """

    # --- Define constants and given values ---
    
    # Reduced Planck constant (J·s)
    hbar = 1.054571817e-34
    
    # Given velocity of the electron (m/s)
    v = 2e8
    
    # Given uncertainty in position (m)
    # Δx = 0.1 nm = 0.1 * 10^-9 m = 1e-10 m
    delta_x = 0.1e-9

    # The answer provided by the LLM to be checked
    llm_answer_option = 'C'

    # --- Perform the calculation ---

    # 1. Calculate the minimum uncertainty in momentum (Δp_x) using Heisenberg's Uncertainty Principle.
    # The principle states Δx * Δp_x ≥ ħ / 2.
    # The minimum uncertainty is therefore Δp_x_min = ħ / (2 * Δx).
    try:
        delta_p_min = hbar / (2 * delta_x)
    except ZeroDivisionError:
        return "Constraint Error: The uncertainty in position (Δx) cannot be zero."

    # 2. Estimate the minimum uncertainty in energy (ΔE).
    # The uncertainty in energy can be related to the uncertainty in momentum via ΔE ≈ v * Δp.
    # This relation holds for both non-relativistic (E = p^2/2m) and relativistic (E^2 = (pc)^2 + (m_0c^2)^2)
    # cases, since dE/dp = v in both regimes.
    delta_E_min = v * delta_p_min

    # --- Verify the answer ---

    # Define the numerical values for the given options
    options = {
        "A": 1e-17,
        "B": 1e-19,
        "C": 1e-16,
        "D": 1e-18
    }

    # Find which option is numerically closest to our calculated result
    calculated_closest_option = min(options.keys(), key=lambda k: abs(options[k] - delta_E_min))

    # Compare the calculated best option with the LLM's answer
    if calculated_closest_option == llm_answer_option:
        return "Correct"
    else:
        # If the LLM's answer is not the closest option, it's incorrect.
        # Provide a detailed reason.
        reason = (
            f"The provided answer '{llm_answer_option}' is incorrect.\n"
            f"The calculation based on physical principles is as follows:\n"
            f"1. The Heisenberg Uncertainty Principle gives the minimum momentum uncertainty: "
            f"Δp_min = ħ / (2 * Δx).\n"
            f"   Δp_min = {hbar:.4e} J·s / (2 * {delta_x:.1e} m) ≈ {delta_p_min:.4e} kg·m/s.\n"
            f"2. The corresponding minimum energy uncertainty is: ΔE_min = v * Δp_min.\n"
            f"   ΔE_min = {v:.1e} m/s * {delta_p_min:.4e} kg·m/s ≈ {delta_E_min:.4e} J.\n"
            f"3. The calculated value {delta_E_min:.4e} J is approximately 1.05 x 10^-16 J.\n"
            f"4. Comparing this result to the options, the calculated value is closest to option C ({options['C']:.0e} J), not the provided option {llm_answer_option} ({options[llm_answer_option]:.0e} J)."
        )
        return reason

# You can run this function to get the verification result.
# For example:
# result = check_correctness()
# print(result)