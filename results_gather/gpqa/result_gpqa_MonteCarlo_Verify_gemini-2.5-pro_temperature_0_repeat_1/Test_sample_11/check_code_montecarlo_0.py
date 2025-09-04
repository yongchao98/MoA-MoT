import math

def check_pion_decay_answer():
    """
    Checks the correctness of the proposed answer for pion decay by verifying
    the conservation of energy and momentum.
    """
    # --- Problem Constants ---
    # Rest masses are given in MeV, so we can treat them as energies (in units of MeV/c^2 * c^2).
    m_pi = 139.6  # Rest energy of Pi(+) in MeV
    m_mu = 105.7  # Rest energy of mu(+) in MeV

    # --- Proposed Answer from LLM (Option B) ---
    # Kinetic energies of the product particles in MeV.
    ke_mu_ans = 4.12
    ke_nu_ans = 29.8

    # Set a relative tolerance for floating-point comparisons.
    # This accounts for potential rounding in the provided answer choices.
    # A 1% tolerance is reasonable for this problem's precision.
    relative_tolerance = 0.01

    # --- 1. Verification of Energy Conservation ---
    # Initial energy is the rest energy of the stationary pion.
    e_initial = m_pi

    # Final energy is the sum of the total energies of the products.
    # Total energy = Rest Energy + Kinetic Energy
    e_mu_total = m_mu + ke_mu_ans
    # For a massless neutrino, rest energy is 0, so total energy is its kinetic energy.
    e_nu_total = ke_nu_ans
    e_final = e_mu_total + e_nu_total

    if not math.isclose(e_initial, e_final, rel_tol=relative_tolerance):
        reason = (
            f"Incorrect: The answer violates the conservation of energy.\n"
            f"Initial Energy (Pion Rest Energy): {e_initial:.2f} MeV\n"
            f"Final Energy (Muon Total Energy + Neutrino Total Energy): {e_final:.2f} MeV\n"
            f"The values should be equal."
        )
        return reason

    # --- 2. Verification of Momentum Conservation ---
    # Since the initial pion is at rest, the final momenta must be equal in magnitude.
    # We calculate the momentum of each particle (in units of MeV/c).

    # Muon momentum: (pc)^2 = E^2 - (mc^2)^2
    # We must check for a non-negative value before taking the square root.
    pc_mu_squared = e_mu_total**2 - m_mu**2
    if pc_mu_squared < 0:
        reason = (
            f"Incorrect: The proposed kinetic energy for the muon is physically impossible.\n"
            f"The total energy squared ({e_mu_total**2:.2f}) is less than the rest energy squared ({m_mu**2:.2f}), "
            f"which would imply an imaginary momentum."
        )
        return reason
    pc_mu = math.sqrt(pc_mu_squared)

    # Neutrino momentum: pc = E (since it's massless)
    pc_nu = e_nu_total

    if not math.isclose(pc_mu, pc_nu, rel_tol=relative_tolerance):
        reason = (
            f"Incorrect: The answer violates the conservation of momentum.\n"
            f"The product particles must have momenta of equal magnitude.\n"
            f"Muon Momentum Magnitude (p*c): {pc_mu:.2f} MeV\n"
            f"Neutrino Momentum Magnitude (p*c): {pc_nu:.2f} MeV\n"
            f"The values should be equal."
        )
        return reason

    # If both checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_pion_decay_answer()
print(result)