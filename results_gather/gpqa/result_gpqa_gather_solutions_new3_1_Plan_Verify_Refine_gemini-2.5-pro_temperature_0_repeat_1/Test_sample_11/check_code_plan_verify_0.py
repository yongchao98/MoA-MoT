import math

def check_correctness():
    """
    Checks the correctness of the answer for the pion decay problem.

    The problem involves a two-body decay of a stationary particle (pion)
    into a massive particle (muon) and a massless particle (neutrino).
    The solution must satisfy two main constraints:
    1. Conservation of Energy: The total kinetic energy of the products must equal
       the difference in rest mass energy between the initial and final particles (the Q-value).
    2. Conservation of Momentum: The product particles must have momenta of equal magnitude.
    """
    # Given values from the question
    m_pi_c2 = 139.6  # Rest mass energy of Pi(+) in MeV
    m_mu_c2 = 105.7  # Rest mass energy of mu(+) in MeV

    # The final answer from the LLM's analysis is D.
    # The options are:
    # A) 7.2 MeV, 32.8 MeV
    # B) 3.52 MeV, 20.8 MeV
    # C) 2.84 MeV, 26.8 MeV
    # D) 4.12 MeV, 29.8 MeV
    
    # From the physics, the lighter particle (neutrino) carries more energy.
    # The heavier particle (muon) has less kinetic energy.
    # So, for option D, we assign:
    ke_mu_proposed = 4.12  # MeV
    ke_nu_proposed = 29.8  # MeV

    # A small tolerance is needed for floating-point comparisons due to rounding in the options.
    tolerance = 0.05

    # --- Constraint 1: Conservation of Energy ---
    # The total kinetic energy released (Q-value) must equal the difference in rest mass energies.
    # Q = KE_mu + KE_nu = m_pi*c^2 - m_mu*c^2
    q_value_expected = m_pi_c2 - m_mu_c2
    q_value_proposed = ke_mu_proposed + ke_nu_proposed

    if not math.isclose(q_value_expected, q_value_proposed, abs_tol=tolerance):
        return (f"Incorrect. The answer does not satisfy the conservation of energy. "
                f"The sum of the kinetic energies should be the Q-value. "
                f"Expected sum (Q-value): {q_value_expected:.2f} MeV. "
                f"Proposed sum from option D: {q_value_proposed:.2f} MeV.")

    # --- Constraint 2: Conservation of Momentum ---
    # The magnitudes of the momenta of the two product particles must be equal (p_mu = p_nu).
    
    # For the massless neutrino, its energy is purely kinetic, and E_nu = p_nu*c.
    # So, (p_nu*c)^2 = (KE_nu)^2
    p_nu_c_sq = ke_nu_proposed**2

    # For the massive muon, we use the relativistic energy-momentum relation: E_mu^2 = (p_mu*c)^2 + (m_mu*c^2)^2
    # where the total energy of the muon is E_mu = KE_mu + m_mu*c^2.
    # So, (p_mu*c)^2 = E_mu^2 - (m_mu*c^2)^2 = (KE_mu + m_mu*c^2)^2 - (m_mu*c^2)^2
    e_mu_proposed = ke_mu_proposed + m_mu_c2
    p_mu_c_sq = e_mu_proposed**2 - m_mu_c2**2

    # Check if the squared momenta are equal. Using squared values avoids sqrt and is numerically stable.
    # A larger tolerance is used for squared values.
    if not math.isclose(p_mu_c_sq, p_nu_c_sq, abs_tol=tolerance * 100):
        p_mu_c = math.sqrt(p_mu_c_sq)
        p_nu_c = math.sqrt(p_nu_c_sq)
        return (f"Incorrect. The answer does not satisfy the conservation of momentum. "
                f"The momentum magnitudes of the muon and neutrino should be equal. "
                f"Calculated muon momentum (p*c): {p_mu_c:.2f} MeV. "
                f"Calculated neutrino momentum (p*c): {p_nu_c:.2f} MeV.")

    # If both constraints are satisfied, the answer is correct.
    return "Correct"

# Execute the check
result = check_correctness()
print(result)