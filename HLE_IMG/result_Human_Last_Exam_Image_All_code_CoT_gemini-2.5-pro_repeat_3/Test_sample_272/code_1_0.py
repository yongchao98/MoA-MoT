import math

def calculate_fine_tuning_factor():
    """
    Calculates the multiplicative factor change in the fine-tuning measure.
    """
    # Given physical constants and scales
    m_h_phys = 125.5  # GeV
    y_t = 0.95
    
    # Scales in GeV
    Lambda1 = 8.0 * 1e3   # 8 TeV
    Lambda2 = 1.1 * 1e6   # 1.1 PeV

    # Calculate the coefficient from the loop correction
    coeff = y_t**2 / (16 * math.pi**2)

    # Calculate the arguments of the square roots for the bare mass calculation
    # m_h_bare^2 = m_h_phys^2 + coeff * Lambda^2
    term_Lambda1 = m_h_phys**2 + coeff * Lambda1**2
    term_Lambda2 = m_h_phys**2 + coeff * Lambda2**2

    # The multiplicative factor is the ratio of the bare masses
    factor = math.sqrt(term_Lambda2) / math.sqrt(term_Lambda1)

    print("This script calculates the change in a fine-tuning measure for the Higgs mass.")
    print("The final equation for the factor is derived from the ratio of the bare Higgs mass at two different energy scales.")
    print("\nGiven values:")
    print(f"Physical Higgs mass (m_h_phys): {m_h_phys} GeV")
    print(f"Top Yukawa coupling (y_t): {y_t}")
    print(f"Initial scale (Lambda_1): {int(Lambda1)} GeV")
    print(f"New scale (Lambda_2): {int(Lambda2)} GeV")
    
    print("\nCalculation of the factor:")
    print(f"Factor = sqrt(m_h_phys^2 + (y_t^2 / (16*pi^2)) * Lambda_2^2) / sqrt(m_h_phys^2 + (y_t^2 / (16*pi^2)) * Lambda_1^2)")
    print(f"Factor = sqrt({m_h_phys:.1f}^2 + ({y_t:.2f}^2 / (16*pi^2)) * {Lambda2:.1e}^2) / sqrt({m_h_phys:.1f}^2 + ({y_t:.2f}^2 / (16*pi^2)) * {Lambda1:.1e}^2)")
    print(f"Factor = sqrt({term_Lambda2:.4e}) / sqrt({term_Lambda1:.4e})")
    print(f"Factor = {math.sqrt(term_Lambda2):.2f} / {math.sqrt(term_Lambda1):.2f}")
    
    print("\nFinal Result:")
    print(f"The multiplicative factor by which the fine-tuning changes is: {factor:.2f}")

calculate_fine_tuning_factor()
<<<134.64>>>