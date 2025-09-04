import math

def check_particle_decay_ke():
    """
    Checks the kinetic energies of product particles in the decay Pi(+) -> mu(+) + nu.

    The function calculates the theoretical kinetic energies based on conservation laws
    and compares them to the values given in the proposed correct answer (Option C).
    """
    # --- Given Constants ---
    # Rest mass energy of the positive pion (Pi+) in MeV
    m_pi_c2 = 139.6
    # Rest mass energy of the positive muon (mu+) in MeV
    m_mu_c2 = 105.7
    # The neutrino (nu) is assumed to be massless.

    # --- Proposed Answer Values ---
    # The final answer selected is C, which corresponds to the values:
    # KE_mu = 4.12 MeV, KE_nu = 29.8 MeV
    # The order doesn't matter for the check, so we can sort them.
    answer_ke_values = sorted([4.12, 29.8])

    # --- Theoretical Calculation ---
    # This calculation is based on the principles of conservation of energy and momentum
    # for a relativistic two-body decay from a stationary particle.

    # 1. Calculate the kinetic energy of the massless neutrino (KE_nu).
    # The formula is derived from the conservation laws:
    # KE_nu = pc = ((m_pi_c2)^2 - (m_mu_c2)^2) / (2 * m_pi_c2)
    try:
        ke_nu = (m_pi_c2**2 - m_mu_c2**2) / (2 * m_pi_c2)

        # 2. Calculate the total kinetic energy released (Q-value).
        # This is the difference between the initial and final rest mass energies.
        q_value = m_pi_c2 - m_mu_c2

        # 3. Calculate the kinetic energy of the massive muon (KE_mu).
        # The total kinetic energy is shared between the two particles.
        ke_mu = q_value - ke_nu

        calculated_ke_values = sorted([ke_mu, ke_nu])

    except Exception as e:
        return f"An error occurred during the theoretical calculation: {e}"

    # --- Verification ---
    # Compare the calculated values with the values from the proposed answer.
    # A small tolerance is used to account for potential rounding in the answer.
    tolerance = 0.01

    # Check if the individual calculated values match the answer's values.
    is_muon_ke_correct = math.isclose(calculated_ke_values[0], answer_ke_values[0], rel_tol=tolerance)
    is_neutrino_ke_correct = math.isclose(calculated_ke_values[1], answer_ke_values[1], rel_tol=tolerance)

    if is_muon_ke_correct and is_neutrino_ke_correct:
        return "Correct"
    else:
        reason = "The answer is incorrect.\n"
        reason += f"Calculated KE for muon: {calculated_ke_values[0]:.3f} MeV, Answer's value: {answer_ke_values[0]} MeV.\n"
        reason += f"Calculated KE for neutrino: {calculated_ke_values[1]:.3f} MeV, Answer's value: {answer_ke_values[1]} MeV.\n"
        
        # Also check the conservation of energy as a constraint.
        answer_q_value = sum(answer_ke_values)
        if not math.isclose(answer_q_value, q_value, rel_tol=tolerance):
            reason += f"Constraint Violated: The sum of kinetic energies in the answer ({answer_q_value:.2f} MeV) does not equal the total energy released (Q-value), which should be {q_value:.2f} MeV."
        
        return reason

# Run the check
result = check_particle_decay_ke()
print(result)