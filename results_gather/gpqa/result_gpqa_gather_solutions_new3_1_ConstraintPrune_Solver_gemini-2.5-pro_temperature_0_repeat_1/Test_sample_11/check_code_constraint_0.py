import math

def check_correctness():
    """
    Checks if the proposed kinetic energies for the pion decay satisfy
    the laws of conservation of energy and momentum.
    """
    # --- Given values from the question ---
    # All energies are in MeV.
    m_pi = 139.6  # Rest mass energy of the pion
    m_mu = 105.7  # Rest mass energy of the muon
    m_nu = 0.0    # Assumed rest mass energy of the neutrino

    # --- Proposed answer from the LLM (Option A) ---
    # The first value corresponds to the massive particle (muon), 
    # the second to the massless one (neutrino).
    ke_mu_answer = 4.12
    ke_nu_answer = 29.8

    # --- Constraint 1: Conservation of Energy ---
    # The total kinetic energy released (Q-value) must equal the difference 
    # in rest mass energy between the initial and final states.
    
    # Theoretical Q-value
    q_value_theoretical = m_pi - (m_mu + m_nu)
    
    # Q-value from the proposed answer
    q_value_answer = ke_mu_answer + ke_nu_answer

    # Check if the Q-values match, allowing for small rounding differences.
    if not math.isclose(q_value_theoretical, q_value_answer, rel_tol=1e-3):
        return (f"Incorrect. The conservation of energy is not satisfied. "
                f"The sum of the kinetic energies in the answer is {q_value_answer:.2f} MeV, "
                f"but the theoretical total kinetic energy (Q-value) released should be "
                f"{q_value_theoretical:.2f} MeV.")

    # --- Constraint 2: Conservation of Momentum ---
    # Since the initial pion is at rest, the final particles must have equal
    # and opposite momenta. We check if their momentum magnitudes are equal.
    # We will compare the square of the momentum, (pc)^2, for each particle.

    # For the massive muon: (pc)^2 = E^2 - (mc^2)^2
    # where E = KE + mc^2
    total_energy_mu = ke_mu_answer + m_mu
    pc_squared_mu = total_energy_mu**2 - m_mu**2
    
    # For the massless neutrino: E = pc
    # Since E = KE for a massless particle, pc = KE.
    pc_squared_nu = ke_nu_answer**2
    
    # Check if the squared momenta match, allowing for small rounding differences.
    if not math.isclose(pc_squared_mu, pc_squared_nu, rel_tol=1e-3):
        return (f"Incorrect. The conservation of momentum is not satisfied. "
                f"The magnitudes of the momenta of the two particles are not equal. "
                f"Calculated (p_mu*c)^2 = {pc_squared_mu:.2f} MeV^2. "
                f"Calculated (p_nu*c)^2 = {pc_squared_nu:.2f} MeV^2.")

    # If both constraints are satisfied, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)