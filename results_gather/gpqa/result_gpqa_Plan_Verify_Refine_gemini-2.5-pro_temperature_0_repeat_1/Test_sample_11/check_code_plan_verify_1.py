import math

def check_particle_decay_ke():
    """
    This function verifies the kinetic energies of the products in the decay of a stationary pion.
    Decay: Pi(+) -> mu(+) + nu
    The function uses the principles of conservation of energy and momentum in relativistic mechanics.
    """

    # --- Given Constants ---
    # Rest energy of the positive pion (Pi+) in MeV
    m_pi_c2 = 139.6
    # Rest energy of the positive muon (mu+) in MeV
    m_mu_c2 = 105.7
    # The neutrino (nu) is assumed to be massless, which is a standard approximation in this context.
    m_nu_c2 = 0

    # --- Answer to be Checked (Option D) ---
    # The option provides a pair of kinetic energies for the product particles.
    # We use a set to be independent of the order.
    expected_ke_values = {4.12, 29.8}  # in MeV

    # --- Physics Calculation ---
    # For a stationary particle decaying into two products, the kinetic energy of product 1 (muon) is given by:
    # KE_1 = ( (M*c^2 - m1*c^2 - m2*c^2)^2 + 2*(M*c^2)*(m2*c^2 - m1*c^2) + (m1*c^2 - m2*c^2)^2 ) / (2 * M*c^2)
    # A simpler formula exists when one product (neutrino, m2) is massless:
    # KE_mu = ( (m_pi*c^2 - m_mu*c^2)^2 ) / ( 2 * m_pi*c^2 )
    
    try:
        # Calculate the kinetic energy of the muon
        calculated_ke_mu = ((m_pi_c2 - m_mu_c2)**2) / (2 * m_pi_c2)

        # The total kinetic energy released (Q-value) is the difference in rest mass energy
        q_value = m_pi_c2 - m_mu_c2 - m_nu_c2

        # The kinetic energy of the neutrino is the remaining energy
        calculated_ke_nu = q_value - calculated_ke_mu
        
        calculated_values = {calculated_ke_mu, calculated_ke_nu}

    except ZeroDivisionError:
        return "Calculation Error: The rest mass of the initial particle (pion) cannot be zero."
    except Exception as e:
        return f"An unexpected error occurred during calculation: {e}"

    # --- Verification ---
    # Compare the calculated set of energies with the expected set from the answer.
    # A tolerance is used to account for rounding in the provided answer.
    # The answer is given to two decimal places, so a tolerance of 0.02 is reasonable.
    tolerance = 0.02
    
    # Unpack the sets for comparison
    expected_val1, expected_val2 = sorted(list(expected_ke_values))
    calc_val1, calc_val2 = sorted(list(calculated_values))

    if math.isclose(calc_val1, expected_val1, abs_tol=tolerance) and \
       math.isclose(calc_val2, expected_val2, abs_tol=tolerance):
        
        # Sanity check: The heavier particle (muon) must have lower kinetic energy.
        if calculated_ke_mu > calculated_ke_nu:
            return (f"Incorrect. Although the energy values are correct, they are assigned incorrectly. "
                    f"The heavier muon should have less kinetic energy. "
                    f"Calculated K_mu={calculated_ke_mu:.3f} MeV, K_nu={calculated_ke_nu:.3f} MeV.")
        
        return "Correct"
    else:
        # If the check fails, provide a detailed reason.
        reason = "Incorrect. The calculated kinetic energies do not match the provided answer.\n"
        reason += f"Provided answer values (sorted): {expected_val1} MeV, {expected_val2} MeV\n"
        reason += f"Calculated values (sorted): {calc_val1:.4f} MeV, {calc_val2:.4f} MeV\n"
        reason += f"The difference between calculated and expected values exceeds the tolerance of {tolerance}."
        return reason

# Execute the check and print the result.
result = check_particle_decay_ke()
# The code will return "Correct" if the logic and values in the LLM's answer are valid.
# Otherwise, it will return a string explaining the error.
# In this case, the LLM's answer is correct.
print(result)