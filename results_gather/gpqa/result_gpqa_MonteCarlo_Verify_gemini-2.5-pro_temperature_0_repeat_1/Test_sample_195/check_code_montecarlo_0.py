import numpy as np

def check_correctness_of_relativistic_vmax():
    """
    This function checks the correctness of the provided answer for the maximum speed
    of a 1D relativistic harmonic oscillator.

    The question is:
    Consider a 1-dimensional relativistic harmonic oscillator with mass m and maximum
    amplitude A obeying Hook's law (F=-kx). What is the maximum speed v_max of the mass?
    The speed of light is c.

    The proposed answer is C:
    v_max = c * sqrt(1 - 1 / (1 + (k*A^2)/(2*m*c^2))^2)
    """

    # --- Step 1: Define the ground truth formula from first principles ---
    # The total energy E is conserved.
    # At maximum displacement (x=A), velocity v=0. The energy is purely potential.
    # E_total = U_max + K(v=0) = (1/2)kA^2 + (gamma(0)-1)mc^2 = (1/2)kA^2
    #
    # At the equilibrium position (x=0), velocity is maximum (v=v_max). The energy is purely kinetic.
    # E_total = U(x=0) + K(v=v_max) = 0 + (gamma_max - 1)mc^2
    #
    # By conservation of energy:
    # (1/2)kA^2 = (gamma_max - 1)mc^2
    #
    # Solve for gamma_max:
    # gamma_max = 1 + (k * A^2) / (2 * m * c^2)
    #
    # Now, solve for v_max from the definition of gamma_max = 1 / sqrt(1 - v_max^2/c^2):
    # v_max = c * sqrt(1 - 1 / gamma_max^2)
    #
    # Substituting the expression for gamma_max gives the ground truth formula.
    def ground_truth_vmax(m, k, A, c):
        try:
            energy_ratio = (k * A**2) / (2 * m * c**2)
            gamma_max = 1 + energy_ratio
            
            # The term inside the sqrt must be non-negative. Since gamma_max >= 1,
            # this condition is always met for physical parameters.
            term_in_sqrt = 1 - 1 / (gamma_max**2)
            
            return c * np.sqrt(term_in_sqrt)
        except (ValueError, OverflowError, ZeroDivisionError):
            return np.nan

    # --- Step 2: Define the formula from the proposed answer (Option C) ---
    def answer_vmax(m, k, A, c):
        try:
            denominator_term = 1 + (k * A**2) / (2 * m * c**2)
            term_in_sqrt = 1 - 1 / (denominator_term**2)
            return c * np.sqrt(term_in_sqrt)
        except (ValueError, OverflowError, ZeroDivisionError):
            return np.nan

    # --- Step 3: Perform numerical verification with random parameters ---
    c_light = 299792458.0
    num_trials = 2000
    for i in range(num_trials):
        # Generate random physical parameters over several orders of magnitude
        m = 10**np.random.uniform(-30, 2)  # Mass in kg (electron to person)
        k = 10**np.random.uniform(-5, 10)   # Spring constant in N/m
        A = 10**np.random.uniform(-15, 1)  # Amplitude in m (atomic to macroscopic)

        truth_val = ground_truth_vmax(m, k, A, c_light)
        answer_val = answer_vmax(m, k, A, c_light)

        if np.isnan(truth_val) or np.isnan(answer_val):
            continue # Skip trials that result in numerical errors from extreme values

        # Compare the results. They should be identical within a very small tolerance.
        if not np.isclose(truth_val, answer_val, rtol=1e-12, atol=1e-15):
            return (f"Incorrect. The answer's formula does not match the first-principles derivation.\n"
                    f"For a trial with parameters m={m:.2e}, k={k:.2e}, A={A:.2e}:\n"
                    f"  - Ground Truth v_max: {truth_val}\n"
                    f"  - Answer's v_max:     {answer_val}")

    # --- Step 4: Check the non-relativistic limit ---
    # In this limit (kA^2 << mc^2), the result should approach the classical v_max = sqrt(k*A^2/m).
    m_nr, k_nr, A_nr = 1.0, 1.0, 1e-7
    v_classical = np.sqrt(k_nr * A_nr**2 / m_nr)
    v_answer_nr = answer_vmax(m_nr, k_nr, A_nr, c_light)
    if not np.isclose(v_classical, v_answer_nr):
        return (f"Incorrect. The answer fails the non-relativistic limit check.\n"
                f"  - Classical v_max: {v_classical}\n"
                f"  - Answer's v_max:  {v_answer_nr}")

    # --- Step 5: Check the physical constraint v_max < c ---
    # The structure of the formula c * sqrt(1 - positive_term) where the sqrt is < 1
    # inherently satisfies this constraint. No explicit check is needed if the derivation is correct.

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_correctness_of_relativistic_vmax()
print(result)