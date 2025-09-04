import math

def check_correctness():
    """
    This function checks the correctness of the provided answer for the particle physics decay problem.
    It calculates the kinetic energies of the product particles based on relativistic conservation laws
    and compares them to the values given in the selected answer option.
    """

    # --- Given Constants ---
    # Rest mass energy of the positive pion (Pi+) in MeV
    m_pi_c2 = 139.6
    # Rest mass energy of the positive muon (mu+) in MeV
    m_mu_c2 = 105.7
    # The neutrino is assumed to be massless.

    # --- Answer to be Checked ---
    # The provided answer is B, which corresponds to the values 4.12 MeV and 29.8 MeV.
    # In a two-body decay from rest, the lighter particle (neutrino) gets more kinetic energy.
    # So, KE_mu = 4.12 MeV and KE_nu = 29.8 MeV.
    expected_ke_mu = 4.12
    expected_ke_nu = 29.8

    # --- Physics Calculation ---
    # The calculation is based on the conservation of energy and momentum in special relativity.
    # For a two-body decay of a particle at rest, the energy of the neutrino (particle 2, assumed massless) is:
    # E_nu = KE_nu = pc = ((m_pi*c^2)^2 - (m_mu*c^2)^2) / (2 * m_pi*c^2)
    
    try:
        # Calculate the kinetic energy of the neutrino
        ke_nu_calc = (m_pi_c2**2 - m_mu_c2**2) / (2 * m_pi_c2)

        # The total kinetic energy released (Q-value) is the difference in rest mass energies
        ke_total = m_pi_c2 - m_mu_c2

        # The kinetic energy of the muon is the total KE minus the neutrino's KE
        ke_mu_calc = ke_total - ke_nu_calc
    except ZeroDivisionError:
        return "Error: Division by zero in the calculation. The mass of the pion cannot be zero."
    except Exception as e:
        return f"An unexpected error occurred during calculation: {e}"

    # --- Verification ---
    # Compare the calculated values with the expected values from the answer.
    # A tolerance is used for floating-point comparisons to account for rounding in the options.
    # A tolerance of 0.01 (or 1%) is reasonable for values rounded to two decimal places.
    
    muon_check = math.isclose(ke_mu_calc, expected_ke_mu, abs_tol=0.01)
    neutrino_check = math.isclose(ke_nu_calc, expected_ke_nu, abs_tol=0.01)

    if muon_check and neutrino_check:
        return "Correct"
    else:
        reason = "Incorrect. The calculated kinetic energies do not match the values in the selected answer.\n"
        reason += f"Calculation based on conservation laws gives:\n"
        reason += f"  - KE of muon (mu+): {ke_mu_calc:.2f} MeV\n"
        reason += f"  - KE of neutrino (nu): {ke_nu_calc:.2f} MeV\n"
        reason += f"The selected answer provides:\n"
        reason += f"  - KE of muon (mu+): {expected_ke_mu} MeV\n"
        reason += f"  - KE of neutrino (nu): {expected_ke_nu} MeV\n"
        if not muon_check:
            reason += "The value for the muon's kinetic energy is incorrect.\n"
        if not neutrino_check:
            reason += "The value for the neutrino's kinetic energy is incorrect.\n"
        return reason

# The final response will be the output of this function.
# print(check_correctness())