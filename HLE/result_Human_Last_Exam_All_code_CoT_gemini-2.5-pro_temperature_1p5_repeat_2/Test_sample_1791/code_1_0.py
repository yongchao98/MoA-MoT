import math

def derive_effective_interaction():
    """
    This function explains and prints the derivation of the effective 
    electron-electron interaction by integrating out phonon modes.
    """

    # Introduction to the physical problem
    print("Deriving the effective electron-electron interaction from electron-phonon coupling.")
    print("--------------------------------------------------------------------------------\n")
    print("The goal is to eliminate the phonon degrees of freedom from the total Hamiltonian")
    print("to find an effective interaction potential between the electrons themselves.")
    print("This is done by performing a Gaussian path integral over the phonon fields.\n")
    
    print("The procedure results in an effective interaction Hamiltonian for a specific mode 'q'")
    print("of the form: H_eff = U(q) * rho_q * rho_-q\n")
    print("Here, rho_q is the electron density operator at wavevector q, and U(q) is the")
    print("effective interaction potential mediated by the exchange of a virtual phonon.\n")

    # Define the symbolic components of the final equation
    # These are represented as strings for printing
    g_sq = "g**2"
    qj_sq = "q_j**2"
    m_val = "m"
    w_q_sq = "w_q**2"
    rho_q = "rho_q"
    rho_mq = "rho_-q"

    # Explain the result for the potential U(q)
    print("In the static limit (zero frequency exchange), the potential U(q) is found to be:\n")

    # Print the potential U(q) formulaically
    numerator = f"- ({g_sq} * {qj_sq})"
    denominator = f"({m_val} * {w_q_sq})"
    bar_length = max(len(numerator), len(denominator))
    bar = '-' * bar_length

    print(f"U(q) =  {numerator}")
    print(f"        {bar}")
    print(f"        {denominator}\n")
    print("The negative sign indicates that this phonon-mediated interaction is attractive.\n")

    # Present the final H_eff with each component clearly shown
    print("Therefore, the complete effective electron-electron interaction Hamiltonian for mode q,")
    print("with each symbolic component explicitly shown, is:\n")
    
    print("H_eff = U(q) * rho_q * rho_-q")
    print("      |")
    print("      +--> Broken down into its components:")
    
    # Print the equation with bracketed components
    # The components are defined above as strings
    print(f"H_eff = [ -({g_sq})*({qj_sq}) / (({m_val})*({w_q_sq})) ] * ({rho_q}) * ({rho_mq})")


if __name__ == '__main__':
    derive_effective_interaction()
    
    # Final answer in the required format
    # The result is the expression for the effective Hamiltonian.
    final_expression = "H_eff = - (g**2 * q_j**2) / (m * w_q**2) * rho_q * rho_-q"
    print(f"\n<<<answer\n{final_expression}\n>>>")
