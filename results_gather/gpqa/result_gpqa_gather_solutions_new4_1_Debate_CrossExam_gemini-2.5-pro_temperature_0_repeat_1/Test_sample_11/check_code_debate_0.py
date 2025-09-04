import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the provided answer for the physics problem.
    It calculates the kinetic energies of the decay products based on relativistic conservation laws
    and compares them to the values given in the selected answer.
    """

    # --- Define constants and answer values ---

    # Given rest mass energies from the question
    m_pi_c2 = 139.6  # MeV
    m_mu_c2 = 105.7  # MeV
    # The neutrino is assumed to be massless
    m_nu_c2 = 0.0    # MeV

    # The selected answer is A, which corresponds to the values 4.12 MeV and 29.8 MeV.
    # The detailed explanation identifies 4.12 MeV as the muon's KE and 29.8 MeV as the neutrino's KE.
    answer_ke_mu = 4.12  # MeV
    answer_ke_nu = 29.8  # MeV

    # --- Perform the calculation from first principles ---

    # The kinetic energy of the massless neutrino (KE_nu) can be derived from
    # conservation of energy and momentum.
    # Formula: KE_nu = [ (m_pi*c^2)^2 - (m_mu*c^2)^2 ] / [ 2 * m_pi*c^2 ]
    try:
        calculated_ke_nu = (m_pi_c2**2 - m_mu_c2**2) / (2 * m_pi_c2)
    except ZeroDivisionError:
        return "Calculation failed: The rest mass of the initial particle (pion) cannot be zero."

    # The total kinetic energy released (Q-value) is the difference in rest mass energies.
    # Q = (m_pi - m_mu - m_nu) * c^2
    q_value = m_pi_c2 - m_mu_c2 - m_nu_c2

    # This total kinetic energy is shared between the two products: Q = KE_mu + KE_nu
    # Therefore, the kinetic energy of the muon is:
    calculated_ke_mu = q_value - calculated_ke_nu

    # --- Compare calculated values with the answer ---

    # A tolerance is needed for comparing floating-point numbers, especially
    # since the answer values are rounded. A tolerance of 0.02 should be
    # sufficient to account for rounding differences.
    tolerance = 0.02

    # Check if the calculated kinetic energy for the muon matches the answer
    muon_ke_matches = abs(calculated_ke_mu - answer_ke_mu) < tolerance

    # Check if the calculated kinetic energy for the neutrino matches the answer
    neutrino_ke_matches = abs(calculated_ke_nu - answer_ke_nu) < tolerance

    # --- Return the result ---

    if muon_ke_matches and neutrino_ke_matches:
        return "Correct"
    else:
        error_message = "Incorrect. The calculated kinetic energies do not match the answer's values.\n"
        error_message += f"Calculation based on conservation laws:\n"
        error_message += f"  - Calculated KE of muon: {calculated_ke_mu:.4f} MeV\n"
        error_message += f"  - Calculated KE of neutrino: {calculated_ke_nu:.4f} MeV\n"
        error_message += f"Answer's values (Option A):\n"
        error_message += f"  - KE of muon: {answer_ke_mu} MeV\n"
        error_message += f"  - KE of neutrino: {answer_ke_nu} MeV\n"
        
        if not muon_ke_matches:
            error_message += f"Reason: The calculated muon KE ({calculated_ke_mu:.4f}) is not within the tolerance ({tolerance}) of the answer's value ({answer_ke_mu}).\n"
        if not neutrino_ke_matches:
            error_message += f"Reason: The calculated neutrino KE ({calculated_ke_nu:.4f}) is not within the tolerance ({tolerance}) of the answer's value ({answer_ke_nu}).\n"
            
        # Check if the values might have been swapped
        muon_swapped_matches = abs(calculated_ke_mu - answer_ke_nu) < tolerance
        neutrino_swapped_matches = abs(calculated_ke_nu - answer_ke_mu) < tolerance
        if muon_swapped_matches and neutrino_swapped_matches:
            error_message += "It appears the kinetic energy values in the answer might be swapped between the muon and the neutrino."

        return error_message

# Execute the check
result = check_correctness_of_answer()
print(result)