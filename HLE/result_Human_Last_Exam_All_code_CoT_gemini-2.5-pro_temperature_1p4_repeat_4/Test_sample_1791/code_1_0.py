import sympy as sp

def derive_effective_interaction():
    """
    This function prints the step-by-step derivation of the effective 
    electron-electron interaction mediated by phonons.
    """

    print("### Derivation of the Effective Electron-Electron Interaction ###\n")
    print("This script outlines the derivation by performing a path integral over the phonon fields.\n")

    # --- Step 1: Define the System ---
    print("--- Step 1: Starting Hamiltonian ---")
    print("The phonon Hamiltonian is given by:")
    print("H_ph = sum_{q, j} w_q * a†_{q,j} * a_{q,j}\n")
    print("The electron-phonon interaction is:")
    print("H_el-ph = g * sum_{k,q,j} (i*q_j / (2*m*w_q)^(1/2)) * n_{q} * (a_{q,j} + a†_{-q,j})")
    print("where 'n_q' is the electron density operator in Fourier space.\n")
    print("Let's define a coupling constant V_{q,j} = g * (i*q_j / (2*m*w_q)^(1/2)).")
    print("The interaction can be rewritten by separating the creation and annihilation operators:")
    print("H_el-ph = sum_{q,j} [ (V_{q,j} * n_q) * a_{q,j} + (V_{-q,j} * n_{-q}) * a†_{q,j} ]\n")

    # --- Step 2: Path Integral Formulation ---
    print("--- Step 2: Path Integral over Phonon Fields ---")
    print("We integrate out the phonon fields (a, a†) using the path integral formalism in imaginary time tau.")
    print("The phonon part of the action is S_ph + S_el-ph.")
    print("S_ph = integral_0^beta d(tau) sum_{q,j} a*_{q,j}(tau) * (d/d(tau) + w_q) * a_{q,j}(tau)")
    print("S_el-ph = integral_0^beta d(tau) H_el-ph(tau)\n")
    print("The electron density n_q(tau) acts as a source for the phonon fields.\n")
    
    # --- Step 3: Gaussian Integration in Frequency Space ---
    print("--- Step 3: Gaussian Integration in Matsubara Frequency Space ---")
    print("The action is quadratic in the phonon fields, so this is a Gaussian integral.")
    print("After Fourier transforming to Matsubara frequencies (w_n), we get an effective action for the electrons.")
    print("The result of integrating out the phonon fields a_{q,j,n} is an effective action S_eff for the source terms (n_q).")
    print("The general result of such an integral gives an effective action of the form:")
    print("S_eff ~ - sum_{q,j,n} J*_{q,j,n} * G_0(q,j,w_n) * J_{q,j,-n}\n")
    print("Where J* and J are the source terms from H_el-ph, and G_0 is the non-interacting phonon Green's function.")
    print("From our H_el-ph, the sources are:")
    print("J*_{q,j,n} = V_{q,j} * n_{q,n}")
    print("J_{q,j,-n} = V_{-q,j} * n_{-q,-n}\n")
    print("The Green's function for the action specified is G_0(q, j, w_n) = 1 / (w_q - i*w_n).\n")
    
    # --- Step 4: Symmetrize and find the Interaction Potential ---
    print("--- Step 4: Symmetrization and Interaction Potential U(q, w_n) ---")
    print("Substituting the sources and Green's function, the action term for a given (q, n) is:")
    print("S_eff(q,n) ~ - sum_j (V_{q,j} * V_{-q,j}) / (w_q - i*w_n) * n_{q,n} * n_{-q,-n}\n")
    print("The full interaction must be real and symmetric. We obtain the full interaction potential U(q, w_n) by symmetrizing the expression, combining the contributions from (q, w_n) and (-q, -w_n).")
    print("U(q, w_n) = -sum_j [ V_{q,j}*V_{-q,j} / (w_q - i*w_n) + V_{-q,j}*V_{q,j} / (w_q + i*w_n) ]")
    print("          = -sum_j V_{q,j}*V_{-q,j} * [ (1 / (w_q - i*w_n)) + (1 / (w_q + i*w_n)) ]")
    print("Combining the fractions gives:")
    print("U(q, w_n) = -sum_j V_{q,j}*V_{-q,j} * (2 * w_q) / (w_q^2 + w_n^2)\n")

    # --- Step 5: Final Result ---
    print("--- Step 5: Final Result for the Effective Interaction ---")
    print("Now, we substitute the definition of V_{q,j}.")
    print("The product V_{q,j} * V_{-q,j} is calculated as:")
    print("V_{q,j} * V_{-q,j} = [g*i*q_j / (2*m*w_q)^(1/2)] * [g*i*(-q_j) / (2*m*w_q)^(1/2)]")
    print("                 = -g^2 * i^2 * q_j^2 / (2 * m * w_q)")
    print("                 = g^2 * q_j^2 / (2 * m * w_q)\n")
    print("Plugging this into the expression for U(q, w_n):")
    print("U(q, w_n) = -sum_j [g^2 * q_j^2 / (2 * m * w_q)] * (2 * w_q) / (w_q^2 + w_n^2)")
    print("The factors of '2' and 'w_q' cancel out.\n")
    
    print("=" * 60)
    print("The resulting effective electron-electron interaction potential for a given q is:")
    print("\n          g^2 * q_j^2")
    print("U(q, w_n) = - sum_j --------------")
    print("                  m*(w_q^2 + w_n^2)\n")
    print("This is the phonon-mediated interaction potential between electrons. The negative sign indicates that the interaction is attractive.")
    print("The full interaction term in the effective action is S_eff = (1/2) * sum_{q,n} U(q, w_n) * n_{q,n} * n_{-q,-n}.")
    print("=" * 60)
    
    # This part satisfies the "output each number" instruction by explicitly mentioning them.
    print("\nNote on the numbers in the final equation:")
    print("The final expression for the interaction potential contains the following numbers:")
    print("- The exponent '2' on the coupling constant 'g'.")
    print("- The exponent '2' on the wavevector component 'q_j'.")
    print("- The exponent '2' on the phonon frequency 'w_q'.")
    print("- The exponent '2' on the Matsubara frequency 'w_n'.")
    print("- An implicit coefficient of '1' in the numerator and denominator.")

if __name__ == "__main__":
    derive_effective_interaction()
    
    final_interaction_potential = "- sum_j (g**2 * q_j**2) / (m * (w_q**2 + w_n**2))"
    # The final answer is the mathematical expression itself.
    # We use a string to represent it as requested by the format.
    print(f"\n<<<answer\n{final_interaction_potential}\n>>>")