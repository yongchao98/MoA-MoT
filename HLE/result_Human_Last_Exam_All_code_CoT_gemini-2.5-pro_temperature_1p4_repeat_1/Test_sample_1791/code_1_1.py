import sympy

def derive_effective_interaction():
    """
    This script symbolically derives the effective electron-electron interaction potential
    mediated by phonons using the provided Hamiltonian terms.
    """
    # Define the symbols used in the derivation
    # g: electron-phonon coupling strength
    # q_j: j-th component of the phonon wavevector q
    # m: electron mass
    # w_q: phonon frequency for wavevector q (omega_q)
    # w_n: bosonic Matsubara frequency (omega_n)
    g, q_j, m = sympy.symbols('g q_j m')
    w_q = sympy.Symbol('w_q', positive=True)
    w_n = sympy.Symbol('w_n')

    # Define the electron-phonon coupling vertex V_{q,j}
    # V_{q,j} = g * i * q_j / sqrt(2 * m * w_q)
    # We need its magnitude squared, |V_{q,j}|^2
    V_qj_sq = (g**2 * q_j**2) / (2 * m * w_q)

    # Define the phonon propagator D(q, omega_n)
    # D(q, w_n) = 2 * w_q / (w_q^2 + w_n^2)
    D_q_wn = (2 * w_q) / (w_q**2 + w_n**2)

    # The effective electron-electron interaction potential U_eff is given by:
    # U_eff = -|V_{q,j}|^2 * D(q, w_n)
    U_eff = -V_qj_sq * D_q_wn

    # Simplify the final expression
    U_eff_simplified = sympy.simplify(U_eff)

    # Print the derivation steps and the final result
    print("Derivation of the Effective Electron-Electron Interaction U_eff(q_j, w_n)")
    print("-" * 65)
    
    # Printing each term of the equation
    print(f"The squared magnitude of the coupling vertex |V_qj|^2 is:")
    print(f"|V_qj|^2 = {V_qj_sq}\n")

    print(f"The phonon propagator D(q, w_n) is:")
    print(f"D(q, w_n) = {D_q_wn}\n")

    print("The effective interaction U_eff = -|V_qj|^2 * D(q, w_n)")
    print(f"U_eff = -({V_qj_sq}) * ({D_q_wn})")
    
    print("\nAfter simplification, the final expression for the effective interaction is:")
    # The final equation with all its parts (interpreted from "output each number in the final equation")
    numerator, denominator = U_eff_simplified.as_numer_denom()
    print(f"U_eff(q_j, w_n) = ({numerator}) / ({denominator})")
    
derive_effective_interaction()