import math

def check_correctness():
    """
    This function checks the correctness of the provided answer to the physics problem.
    It recalculates the mean decay distance from the given parameters and compares it
    to the value in the selected option C.
    """

    # --- Problem Parameters from the question ---
    # Production energy of meson X
    E_X = 8.0  # GeV
    # Mass of meson X
    m_X = 1.2  # GeV/c^2 (or just GeV in natural units where c=1)
    # Decay width of meson X
    Gamma_X_MeV = 320.0  # MeV

    # --- Physical Constants ---
    # Reduced Planck constant * speed of light (hbar*c)
    # We use the value in GeV*m for unit consistency.
    # hbar*c = 197.3269804 MeV*fm = 0.1973269804 GeV*fm
    # 1 fm = 1e-15 m
    hbar_c_GeV_m = 0.1973269804 * 1e-15  # in GeV*m

    # --- Answer to be checked ---
    # The LLM's answer is C, which corresponds to 4.0655 * 10^-15 m
    llm_answer_option = "C"
    llm_answer_value = 4.0655e-15  # meters

    # --- Calculation ---

    # Step 1: Convert all units to be consistent.
    # Convert Gamma_X from MeV to GeV.
    Gamma_X_GeV = Gamma_X_MeV / 1000.0

    # Step 2: Calculate the proper decay length (c*tau_0).
    # This is the decay distance in the particle's rest frame.
    # Formula: c*tau_0 = hbar*c / Gamma_X
    # Constraint: Gamma_X must be positive for a meaningful lifetime.
    if Gamma_X_GeV <= 0:
        return "Constraint not satisfied: The decay width Gamma_X must be a positive value."
    c_tau0 = hbar_c_GeV_m / Gamma_X_GeV

    # Step 3: Calculate the Lorentz factor (gamma).
    # This factor describes relativistic time dilation.
    # Formula: gamma = E_X / m_X
    # Constraint: For a real, moving particle, total energy must exceed rest mass energy.
    if E_X <= m_X:
        return f"Constraint not satisfied: The total energy E_X ({E_X} GeV) must be greater than the rest mass m_X ({m_X} GeV)."
    gamma = E_X / m_X

    # Step 4: Calculate the relativistic momentum factor (beta*gamma).
    # This is needed to find the distance traveled in the lab frame.
    # Formula: beta*gamma = sqrt(gamma^2 - 1)
    beta_gamma = math.sqrt(gamma**2 - 1)

    # Step 5: Calculate the mean decay distance (L) in the lab frame.
    # Formula: L = (beta*gamma) * (c*tau_0)
    calculated_L = beta_gamma * c_tau0

    # --- Verification ---
    # Check if the calculated value is close to the provided answer's value.
    # A small relative tolerance (e.g., 1%) is used to account for potential rounding
    # in the problem's options or constants.
    relative_tolerance = 0.01

    if math.isclose(calculated_L, llm_answer_value, rel_tol=relative_tolerance):
        return "Correct"
    else:
        # If the check fails, provide a detailed reason.
        options = {
            "A": 5.0223e-15,
            "B": 4.0655e-16,
            "C": 4.0655e-15,
            "D": 5.0223e-16,
        }
        # Find the option that is numerically closest to the calculation.
        closest_option = min(options, key=lambda k: abs(options[k] - calculated_L))
        
        reason = (
            f"The answer is incorrect. The selected option is {llm_answer_option} ({llm_answer_value:.4e} m), "
            f"but the calculated mean decay distance is {calculated_L:.4e} m.\n"
            f"The calculated value is closest to option {closest_option} ({options[closest_option]:.4e} m).\n"
            f"Calculation details:\n"
            f"  - Proper decay length (c*tau_0) = hbar*c / Gamma_X = {c_tau0:.4e} m\n"
            f"  - Lorentz factor (gamma) = E_X / m_X = {gamma:.4f}\n"
            f"  - Relativistic factor (beta*gamma) = sqrt(gamma^2 - 1) = {beta_gamma:.4f}\n"
            f"  - Mean decay distance (L) = beta*gamma * c*tau_0 = {calculated_L:.4e} m"
        )
        return reason

# Execute the function and print the result
result = check_correctness()
print(result)