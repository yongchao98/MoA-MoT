import math

def check_answer():
    """
    Checks the correctness of the provided answer for the pion decay problem.
    
    The function first calculates the theoretical kinetic energies of the products
    based on the principles of conservation of energy and momentum. It then
    compares the provided answer against these calculated values.
    """
    
    # --- Given Constants (in MeV) ---
    # Rest energy is equivalent to m*c^2
    m_pi_c2 = 139.6  # Rest energy of Pi(+)
    m_mu_c2 = 105.7  # Rest energy of mu(+)
    # The neutrino is assumed to be massless for this calculation
    m_nu_c2 = 0.0

    # --- The Answer to be Checked (Option A) ---
    # KE_mu, KE_nu in MeV
    answer_ke_mu = 4.12
    answer_ke_nu = 29.8

    # --- Theoretical Calculation ---
    # From conservation of energy and momentum, we can derive the exact KE for the muon:
    # KE_mu = ( (m_pi_c2 - m_mu_c2)^2 ) / (2 * m_pi_c2)
    # This is a standard result for a two-body decay from rest.
    
    try:
        # Calculate the theoretical kinetic energy of the muon
        numerator = (m_pi_c2 - m_mu_c2)**2
        denominator = 2 * m_pi_c2
        theoretical_ke_mu = numerator / denominator

        # Calculate the theoretical kinetic energy of the neutrino using energy conservation
        # Total KE (Q-value) = m_pi_c2 - m_mu_c2
        q_value = m_pi_c2 - m_mu_c2
        theoretical_ke_nu = q_value - theoretical_ke_mu
    except ZeroDivisionError:
        return "Calculation error: The rest mass of the pion cannot be zero."

    # --- Verification Step ---
    # We will check if the provided answer values are close to the theoretical values.
    # A tolerance is used to account for potential rounding in the question's options.
    # A tolerance of 0.1 MeV is reasonable given the precision of the input values.
    tolerance = 0.1

    # Check the muon's kinetic energy
    if not math.isclose(answer_ke_mu, theoretical_ke_mu, abs_tol=tolerance):
        reason = (
            f"Incorrect: The kinetic energy of the muon (mu+) is incorrect.\n"
            f"Based on conservation laws, the calculated KE_mu should be {theoretical_ke_mu:.2f} MeV.\n"
            f"The provided answer gives KE_mu = {answer_ke_mu} MeV."
        )
        return reason

    # Check the neutrino's kinetic energy
    if not math.isclose(answer_ke_nu, theoretical_ke_nu, abs_tol=tolerance):
        reason = (
            f"Incorrect: The kinetic energy of the neutrino (nu) is incorrect.\n"
            f"Based on conservation laws, the calculated KE_nu should be {theoretical_ke_nu:.2f} MeV.\n"
            f"The provided answer gives KE_nu = {answer_ke_nu} MeV."
        )
        return reason
        
    # If both values match the theoretical calculations, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)