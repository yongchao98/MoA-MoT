import math

def check_correctness():
    """
    This function checks the correctness of the given answer to the physics problem.
    
    Problem: Find the kinetic energies (KE) of the product particles in the decay
    Pi(+) -> mu(+) + nu, where the initial Pi(+) is stationary.
    
    Given masses:
    - Rest mass of Pi(+) (m_pi): 139.6 MeV/c^2
    - Rest mass of mu(+) (m_mu): 105.7 MeV/c^2
    - The neutrino (nu) is assumed to be massless.
    
    Provided Answer (Option B):
    - KE_mu = 4.12 MeV
    - KE_nu = 29.8 MeV
    """
    
    # Given values from the question
    m_pi = 139.6  # MeV/c^2
    m_mu = 105.7  # MeV/c^2
    
    # The answer to be checked (Option B)
    # The problem asks for KE of product particles, mu(+) and nu.
    # Based on the calculation, the smaller KE corresponds to the heavier particle (muon).
    answer_ke_mu = 4.12  # MeV
    answer_ke_nu = 29.8   # MeV

    # --- Physics Calculation ---
    # This is a two-body decay of a stationary particle. We use the principles of
    # conservation of energy and conservation of momentum.
    # Let's use units where the speed of light c = 1.
    
    # From conservation of energy and momentum, the magnitude of the momentum 'p'
    # for each of the two daughter particles can be derived as:
    # p = (m_pi^2 - m_mu^2) / (2 * m_pi)
    
    # First, check if the decay is energetically possible.
    if m_pi <= m_mu:
        return "Incorrect. The decay is not possible as the parent particle's mass is not greater than the daughter particle's mass."

    # Calculate the momentum of the products (in units of MeV/c)
    p_magnitude = (m_pi**2 - m_mu**2) / (2 * m_pi)

    # The kinetic energy of the massless neutrino is equal to its total energy, which is p*c (or just 'p' in our units).
    calculated_ke_nu = p_magnitude

    # The kinetic energy of the massive muon is its total energy minus its rest mass energy.
    # Total energy E_mu = sqrt(p^2 + m_mu^2)
    # KE_mu = E_mu - m_mu
    total_energy_mu = math.sqrt(p_magnitude**2 + m_mu**2)
    calculated_ke_mu = total_energy_mu - m_mu

    # --- Verification ---
    # We compare the calculated values with the provided answer.
    # A small tolerance is used to account for potential rounding in the problem's options.
    # A tolerance of 0.02 MeV is appropriate here, as it's less than 0.1% of the larger value
    # and accommodates the rounding difference.
    tolerance = 0.02

    # Check if the calculated KE for the muon is close to the answer's KE for the muon.
    is_mu_correct = math.isclose(calculated_ke_mu, answer_ke_mu, abs_tol=tolerance)
    
    # Check if the calculated KE for the neutrino is close to the answer's KE for the neutrino.
    is_nu_correct = math.isclose(calculated_ke_nu, answer_ke_nu, abs_tol=tolerance)

    if is_mu_correct and is_nu_correct:
        return "Correct"
    else:
        error_message = "Incorrect. The calculated kinetic energies do not match the values from the provided answer within the tolerance.\n"
        error_message += f"Provided Answer (KE_mu, KE_nu): ({answer_ke_mu} MeV, {answer_ke_nu} MeV)\n"
        error_message += f"Calculated Values (KE_mu, KE_nu): ({calculated_ke_mu:.4f} MeV, {calculated_ke_nu:.4f} MeV)\n"
        
        if not is_mu_correct:
            error_message += f"Reason: The calculated muon KE ({calculated_ke_mu:.4f} MeV) is not close enough to the provided value ({answer_ke_mu} MeV).\n"
        if not is_nu_correct:
            error_message += f"Reason: The calculated neutrino KE ({calculated_ke_nu:.4f} MeV) is not close enough to the provided value ({answer_ke_nu} MeV).\n"
            
        return error_message

# The final answer is the output of the check_correctness function.
result = check_correctness()
print(result)