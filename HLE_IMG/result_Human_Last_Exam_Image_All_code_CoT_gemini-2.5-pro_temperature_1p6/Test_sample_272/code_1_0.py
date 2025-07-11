import math

def calculate_fine_tuning_factor():
    """
    Calculates the multiplicative factor change in a fine-tuning measure.
    """
    # Given physical constants
    m_h_phys = 125.5  # GeV
    m_t_phys = 173.2  # GeV
    y_t = 0.95        # Top Yukawa coupling

    # Regularization scales
    # Convert all units to GeV for consistency
    Lambda_1_TeV = 8.0
    Lambda_1_GeV = Lambda_1_TeV * 1e3 # Old scale in GeV
    
    Lambda_2_PeV = 1.1
    Lambda_2_GeV = Lambda_2_PeV * 1e6 # New scale in GeV

    # Calculate the coefficient from the loop calculation
    coeff = y_t**2 / (16 * math.pi**2)
    m_h_phys_sq = m_h_phys**2

    # Calculate the bare mass at the old scale (Lambda_1)
    m_h_bare_sq_1 = m_h_phys_sq + coeff * Lambda_1_GeV**2
    m_h_bare_1 = math.sqrt(m_h_bare_sq_1)

    # Calculate the bare mass at the new scale (Lambda_2)
    m_h_bare_sq_2 = m_h_phys_sq + coeff * Lambda_2_GeV**2
    m_h_bare_2 = math.sqrt(m_h_bare_sq_2)

    # The multiplicative factor is the ratio of the deltas, which simplifies
    # to the ratio of the bare masses.
    factor = m_h_bare_2 / m_h_bare_1

    # Print the final equation with the calculated numbers
    print(f"The analysis shows the multiplicative factor is the ratio of the bare Higgs masses at the two scales.")
    print("Factor = m_h_bare(Lambda_2) / m_h_bare(Lambda_1)")
    print(f"Factor = {m_h_bare_2:.2f} GeV / {m_h_bare_1:.2f} GeV")
    print(f"The multiplicative factor is {factor:.2f}")

calculate_fine_tuning_factor()