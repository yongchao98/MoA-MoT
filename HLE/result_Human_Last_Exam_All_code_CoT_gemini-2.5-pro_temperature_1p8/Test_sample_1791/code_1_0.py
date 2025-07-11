import sympy as sp

def derive_phonon_mediated_interaction():
    """
    Derives and prints the effective electron-electron interaction potential
    mediated by phonons using symbolic mathematics.
    """
    
    # 1. Define symbolic variables for the parameters in the final effective potential.
    # These symbols correspond to the variables in the initial problem.
    g = sp.Symbol('g', real=True, positive=True, doc="The fundamental electron-phonon coupling constant.")
    m = sp.Symbol('m', real=True, positive=True, doc="The mass parameter from the Hamiltonian (related to ionic mass).")
    w_q = sp.Symbol('omega_q', real=True, positive=True, doc="The frequency of a phonon with momentum q.")
    omega = sp.Symbol('omega', real=True, doc="The energy transferred between the interacting electrons.")
    
    # The term sum_j(q_j^2) becomes the magnitude squared of the vector q.
    # We represent this with a single symbol for clarity in the final expression.
    q_mag_sq = sp.Symbol('q^2', real=True, positive=True, doc="The square of the magnitude of the momentum transfer vector q.")

    # 2. Construct the effective potential based on path integral formalism.
    # The derivation (as outlined in the thinking process) yields an effective potential
    # V_eff(q, omega) = - |M_eff|^2 * D(q, omega)
    # where the effective matrix element squared |M_eff|^2 is g^2*|q|^2 / (2*m*w_q)
    # and the phonon propagator D(q, omega) is 2*w_q / (w_q^2 - omega^2).
    # The product simplifies the expression.

    # Numerator of the final expression
    numerator = g**2 * q_mag_sq
    
    # Denominator of the final expression
    denominator = m * (w_q**2 - omega**2)
    
    # The full expression for the effective potential V_eff
    V_eff = -numerator / denominator

    # Create a sympy Equation object for nice printing
    final_equation = sp.Eq(sp.Symbol('V_eff(q, omega)'), V_eff)

    # 3. Print the results in a clear format.
    print("The effective electron-electron interaction potential, derived by integrating out the phonon fields, is:")
    print("-" * 80)
    
    # Use sympy's pretty print to display the final equation
    sp.pprint(final_equation, use_unicode=False)
    
    print("-" * 80)
    print("In the equation above, each variable represents:")
    
    # Print each component of the formula as requested
    print(f"\nV_eff(q, omega) is the effective potential for momentum transfer q and energy transfer omega.")
    print(f"\n{g.name}: {g.doc}")
    print(f"{q_mag_sq.name}: {q_mag_sq.doc}")
    print(f"{m.name}: {m.doc}")
    print(f"{w_q.name}: {w_q.doc}")
    print(f"{omega.name}: {omega.doc}")
    
    print("\nThis result shows that the interaction is attractive (V_eff < 0) when omega < omega_q, which is the basis for conventional superconductivity.")

if __name__ == '__main__':
    derive_phonon_mediated_interaction()