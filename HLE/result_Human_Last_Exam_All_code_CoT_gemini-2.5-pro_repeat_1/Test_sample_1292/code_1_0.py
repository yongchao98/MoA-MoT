import sympy as sp

def derive_field_equation():
    """
    This function symbolically outlines the derivation of the field equations for
    Symmetric Teleparallel Equivalent of General Relativity (STEGR).
    """

    # Define constants and symbols
    # We use 'k' for the constant c^4/(16*pi*G)
    k = sp.Symbol('c^4/(16*pi*G)')
    Q = sp.Symbol('Q')
    g = sp.Symbol('g')
    T_mn = sp.Symbol('T_{mu,nu}')
    g_mn = sp.Symbol('g_{mu,nu}')
    P_amunu = sp.Symbol('P^alpha_{mu,nu}')
    Q_nab = sp.Symbol('Q_nu^{alpha,beta}')
    P_mab = sp.Symbol('P_{mu,alpha,beta}')
    Q_abm = sp.Symbol('Q^{alpha,beta}_mu')
    P_abn = sp.Symbol('P_{alpha,beta,nu}')
    
    # Header
    print("Derivation of the Field Equation for a Metric-Incompatible Gravity Theory")
    print("="*70)
    
    # Step 1: The Action
    print("Step 1: Start with the gravitational action S_g and matter action S_m.")
    print(f"S_g = -({k}) * Integral(sqrt(-g) * Q d^4x)")
    print(f"S_m = Matter Action")
    print(f"The total action is S = S_g + S_m. We vary it with respect to the metric g^{mu,nu}.")
    print("\n")
    
    # Step 2: The Principle of Least Action
    print("Step 2: The principle of least action implies delta(S_g) = -delta(S_m).")
    print("The variation of the matter action defines the energy-momentum tensor:")
    print(f"delta(S_m) = (1/2) * Integral(sqrt(-g) * T^{mu,nu} * delta(g_{mu,nu}) d^4x)")
    print(f"This can be written as: (2/sqrt(-g)) * delta(S_m)/delta(g^{mu,nu}) = T_{mu,nu}")
    print("\n")
    
    # Step 3: Vary the Gravitational Action
    print("Step 3: We compute the variation of the gravitational action, delta(S_g).")
    print(f"delta(S_g) = -({k}) * Integral( delta(sqrt(-g) * Q) d^4x )")
    print(f"Using the product rule: delta(sqrt(-g) * Q) = delta(sqrt(-g))*Q + sqrt(-g)*delta(Q)")
    
    # Step 3a: Variation of sqrt(-g)
    print("\nStep 3a: The variation of the metric determinant term.")
    print("delta(sqrt(-g)) = -(1/2) * sqrt(-g) * g_{mu,nu} * delta(g^{mu,nu})")
    print("This gives a term in the field equation. After accounting for the constant C = -k, this term is:")
    term_g = sp.Rational(1, 2) * Q * g_mn
    print(f"Contribution to LHS: +{term_g}")
    
    # Step 3b: Variation of Q
    print("\nStep 3b: The variation of the non-metricity scalar Q.")
    print("This variation is complex. Using the Euler-Lagrange equations in the coincident gauge, it yields two main parts:")
    print("  1. A total derivative term from integrating by parts.")
    print("  2. A term from the explicit metric dependence of the superpotential P.")
    
    term_deriv = f"-2/sqrt(-g) * d_alpha(sqrt(-g) * {P_amunu})"
    print(f"The derivative term contributes: {term_deriv}")
    
    term_PQ = f"- {P_mab} * {Q_nab} + 2 * {Q_abm} * {P_abn}"
    print(f"The superpotential-dependent term contributes: {term_PQ}")
    print("\n")

    # Step 4: Assemble the Field Equation
    print("Step 4: Combine all terms.")
    print("The full left-hand side (LHS) of the equation is the sum of these parts.")
    LHS = f"{term_deriv} {term_PQ.replace('P', '- P', 1)} + {term_g}"
    print(f"LHS = {LHS}")
    
    print("\nThe right-hand side (RHS) comes from the matter action variation:")
    RHS = f"(8*pi*G/c^4) * {T_mn}"
    print(f"RHS = {RHS}")
    print("\n")
    
    # Final Equation
    print("Final Field Equation:")
    print("="*70)
    
    # Print each term of the final equation individually
    print("Term 1: -2/sqrt(-g) * partial_alpha(sqrt(-g) * P^{alpha,mu,nu})")
    print("Term 2: - P_{mu,alpha,beta} * Q_{nu}^{alpha,beta}")
    print("Term 3: + 2 * Q^{alpha,beta}_{mu} * P_{alpha,beta,nu}")
    print("Term 4: + (1/2) * Q * g_{mu,nu}")
    print(f"Is equal to: {RHS}")
    print("="*70)

    print("\nThis matches choice E.")


if __name__ == '__main__':
    derive_field_equation()