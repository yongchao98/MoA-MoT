import math

def calculate_cross_section(G_F, m_mu, m_e, s):
    """
    Calculates the total cross-section for e- + nu_e_bar -> mu- + nu_mu_bar.

    Args:
        G_F (float): The Fermi constant.
        m_mu (float): The muon mass.
        m_e (float): The electron mass.
        s (float): The center-of-mass energy squared.
    """
    
    # The formula for the total cross-section in the low-energy limit is:
    # sigma = (G_F^2 * (s - m_mu^2)^2) / (2 * pi * s)
    # This is valid only if the process is kinematically allowed.

    # Step 1: Print the formula with the given values.
    # The prompt requires printing the equation with numbers.
    print("Formula for the total cross-section (σ):")
    print(f"σ = (G_F² * (s - m_μ²)²) / (2 * π * s)")
    print("Substituting the given values:")
    # We show the equation with the raw numbers plugged in.
    print(f"σ = ({G_F}² * ({s} - {m_mu}²)²) / (2 * π * {s})")
    
    # Step 2: Check if the reaction is kinematically possible.
    # The center-of-mass energy (sqrt(s)) must be greater than the sum of the rest masses
    # of the final state particles (m_mu + m_numu). Since m_numu is 0, sqrt(s) >= m_mu.
    # A more complete threshold also considers the initial state: sqrt(s) must be >= m_e,
    # and to produce the final state, sqrt(s) must also be >= m_mu.
    # The true physical threshold for the creation of new particles is that
    # the CoM energy must be at least the sum of their masses.
    # In this case, sqrt(s) must be at least m_mu (since m_nu_mu is 0).
    # However, for a realistic process starting from e- and nu_e_bar, the total mass is m_e.
    # The total energy must be >= sum of final particle masses.
    # We will use the most restrictive physical threshold: sqrt(s) >= m_e + m_mu, which would be 
    # more applicable to a process where both are created from pure energy.
    # A simpler threshold sqrt(s) > m_mu ensures the final state is possible. Let's use that.
    
    com_energy = math.sqrt(s)
    threshold_energy = m_mu # Simplified threshold based on final state particle.

    if com_energy < threshold_energy:
        sigma = 0.0
        print(f"\nCenter-of-mass energy sqrt(s) = {com_energy:.3f} is less than the muon mass m_μ = {m_mu:.3f}.")
        print("The process is kinematically forbidden.")
    else:
        # Step 3: Calculate the cross-section if kinematically allowed.
        pi = math.pi
        numerator = (G_F**2 * (s - m_mu**2)**2)
        denominator = (2 * pi * s)
        sigma = numerator / denominator

    # Step 4: Print the final result.
    print("\nFinal Result:")
    print(f"σ = {sigma}")

# Given values from the problem
G_F_val = 1.0
m_mu_val = 1.0
m_e_val = 1.0
s_val = 2.0

calculate_cross_section(G_F_val, m_mu_val, m_e_val, s_val)

# The result is based on the calculation performed in the python code.
# The code determines that the process is kinematically forbidden under the given values
# because sqrt(s=2) = 1.414 is NOT greater than m_mu=1. Oops, my check was wrong.
# Let's re-evaluate. sqrt(2) is indeed > 1. Let's adjust the threshold explanation and code.
# The correct physical threshold for particle production requires that the center-of-mass energy
# must be at least the sum of the rest masses of the final-state particles.
# For mu- and nu_mu_bar, this is m_mu + 0 = m_mu. So sqrt(s) >= m_mu.
# For the given numbers: sqrt(2) = 1.414 which is > 1.0. So the process is allowed.
# Let me regenerate the code block.
