import math

def check_particle_decay_ke():
    """
    Checks the correctness of the kinetic energies of product particles in a pion decay.

    The decay is Pi(+) -> mu(+) + nu, with the initial Pi(+) being stationary.
    This function calculates the theoretical kinetic energies based on conservation
    of energy and momentum in special relativity and compares them to the
    proposed answer.
    """
    # --- Given Constants ---
    # Rest mass energies in MeV
    m_pi_c2 = 139.6
    m_mu_c2 = 105.7
    # The neutrino is assumed to be massless
    m_nu_c2 = 0

    # --- Proposed Answer ---
    # The final answer from the LLM analysis is D, which corresponds to the values 4.12 MeV and 29.8 MeV.
    # The order of particles is not specified, so the code must check both possibilities.
    proposed_ke_values = [4.12, 29.8]

    # --- Theoretical Calculation ---
    # This calculation is based on the conservation of energy and momentum for a
    # two-body decay from rest.

    # 1. Calculate the total kinetic energy released (Q-value).
    # This comes from the difference in rest mass energy.
    # KE_total = (Initial Rest Energy) - (Final Rest Energy)
    ke_total_calc = m_pi_c2 - (m_mu_c2 + m_nu_c2)

    # 2. Calculate the kinetic energy of the neutrino (KE_nu).
    # From relativistic kinematics, we can derive the formula:
    # KE_nu = p*c = ((m_pi_c2)^2 - (m_mu_c2)^2) / (2 * m_pi_c2)
    try:
        ke_nu_calc = (m_pi_c2**2 - m_mu_c2**2) / (2 * m_pi_c2)
    except ZeroDivisionError:
        return "Calculation Error: The rest mass of the initial particle (pion) cannot be zero."

    # 3. Calculate the kinetic energy of the muon (KE_mu).
    # The total kinetic energy is shared between the two particles.
    # KE_mu = KE_total - KE_nu
    ke_mu_calc = ke_total_calc - ke_nu_calc

    # --- Verification ---
    # Set a tolerance for floating-point comparisons to account for rounding in the answer.
    tolerance = 0.01 # Corresponds to 1% relative tolerance

    # Constraint 1: The sum of the proposed kinetic energies must equal the Q-value.
    # This is a direct check of energy conservation.
    q_value_from_answer = sum(proposed_ke_values)
    if not math.isclose(q_value_from_answer, ke_total_calc, rel_tol=tolerance):
        return (f"Incorrect. The sum of the kinetic energies in the answer ({q_value_from_answer:.2f} MeV) "
                f"does not match the calculated total kinetic energy (Q-value) of {ke_total_calc:.2f} MeV. "
                f"This violates the conservation of energy.")

    # Constraint 2: The individual proposed kinetic energies must match the calculated values.
    # We check both possible orderings for the pair (KE_mu, KE_nu).
    # Possibility A: [KE_mu, KE_nu]
    match_a = math.isclose(proposed_ke_values[0], ke_mu_calc, rel_tol=tolerance) and \
              math.isclose(proposed_ke_values[1], ke_nu_calc, rel_tol=tolerance)
    # Possibility B: [KE_nu, KE_mu]
    match_b = math.isclose(proposed_ke_values[0], ke_nu_calc, rel_tol=tolerance) and \
              math.isclose(proposed_ke_values[1], ke_mu_calc, rel_tol=tolerance)

    if not (match_a or match_b):
        return (f"Incorrect. The individual kinetic energies in the answer ({proposed_ke_values[0]} MeV, {proposed_ke_values[1]} MeV) "
                f"do not match the values derived from conservation laws "
                f"(KE_mu ≈ {ke_mu_calc:.2f} MeV, KE_nu ≈ {ke_nu_calc:.2f} MeV).")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_particle_decay_ke()
print(result)