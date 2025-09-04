import math

def check_correctness():
    """
    This function checks the correctness of the given answer for the particle physics problem.

    The question is:
    Consider this annihilation process
    p + p_bar -> 2A+ + 2A-
    The antiproton is slowly moving and m_A*c^2 = 300 MeV.
    What is the velocity of particle A?
    A) 0.77c
    B) 0.96c
    C) 0.91c
    D) 0.86c

    The provided answer is A, which corresponds to a velocity of 0.77c.
    This code will verify if 0.77c is the correct velocity based on physics principles.
    """

    # --- Define constants and given values ---

    # Rest mass energy of a proton (and antiproton) in MeV.
    # Using a standard high-precision value (CODATA 2018).
    m_p_c2 = 938.27208816  # MeV

    # Rest mass energy of particle A in MeV, from the problem statement.
    m_A_c2 = 300.0  # MeV

    # The velocity from the given answer (A) 0.77c.
    # We work with beta = v/c.
    given_beta = 0.77

    # --- Apply Physics Principles ---

    # 1. Assumption from the problem statement:
    # "The antiproton is slowly moving" implies its kinetic energy is negligible.
    # We assume the initial state (proton and antiproton) is effectively at rest.
    # Therefore, the total initial momentum is zero.

    # 2. Conservation of Energy:
    # The total initial energy is the sum of the rest energies of the proton and antiproton.
    E_initial = 2 * m_p_c2

    # The total final energy is the sum of the total energies of the four resulting A particles.
    # By conservation of momentum (initial momentum is zero), the final particles must
    # have a total momentum of zero. By symmetry, they will all have the same speed and energy.
    # E_final = 4 * E_A, where E_A is the total energy of a single particle A.
    #
    # Applying conservation: E_initial = E_final
    # 2 * m_p_c2 = 4 * E_A
    # So, the energy of a single particle A is:
    n_final_particles = 4
    E_A = E_initial / n_final_particles

    # 3. Relativistic Energy Calculation:
    # The total energy of a relativistic particle is E = gamma * m * c^2,
    # where gamma is the Lorentz factor. We can solve for gamma.
    # gamma = E_A / m_A_c2

    # Sanity check: For a reaction to be possible, the energy per particle must be
    # at least its rest mass energy (i.e., gamma >= 1).
    if E_A < m_A_c2:
        return (f"Incorrect: The reaction is physically impossible. The calculated energy per "
                f"particle A ({E_A:.2f} MeV) is less than its rest mass energy ({m_A_c2:.1f} MeV).")

    calculated_gamma = E_A / m_A_c2

    # 4. Solve for Velocity (beta = v/c):
    # The Lorentz factor is gamma = 1 / sqrt(1 - beta^2).
    # Rearranging for beta gives: beta = sqrt(1 - 1/gamma^2).
    try:
        calculated_beta = math.sqrt(1 - 1 / (calculated_gamma**2))
    except ValueError:
        # This case should not be reached due to the sanity check above.
        return "Error: Calculation failed, gamma was less than 1."

    # --- Verification Step ---

    # Compare the calculated beta with the beta from the given answer.
    # A tolerance is used to account for rounding in the problem's options.
    # A tolerance of 0.01 is appropriate for a two-decimal-place answer.
    tolerance = 0.01

    if abs(calculated_beta - given_beta) < tolerance:
        return "Correct"
    else:
        # If the answer is not correct, provide the calculated result for comparison.
        reason = (
            f"Incorrect: The velocity from the given answer ({given_beta}c) does not match the calculated velocity.\n"
            f"Based on conservation of energy:\n"
            f"1. Initial Energy (at rest): E_i = 2 * m_p*c^2 = {E_initial:.2f} MeV.\n"
            f"2. Energy per final particle A: E_A = E_i / 4 = {E_A:.2f} MeV.\n"
            f"3. Lorentz Factor for A: gamma = E_A / m_A*c^2 = {E_A:.2f} / {m_A_c2:.1f} = {calculated_gamma:.4f}.\n"
            f"4. Calculated Velocity: v/c = sqrt(1 - 1/gamma^2) = {calculated_beta:.4f}.\n"
            f"The calculated velocity is approximately {calculated_beta:.2f}c, which does not match the given answer of {given_beta}c."
        )
        return reason

# To use this checker, you would run the function and inspect its return value.
# For example:
# result = check_correctness()
# print(result)