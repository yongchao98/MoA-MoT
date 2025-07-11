import sympy

def solve_sugra_parameters():
    """
    This script calculates the values of alpha^2 and beta based on the
    provided supergravity problem.
    """
    # Define symbols for the constants
    R = sympy.Symbol('R')
    kappa = sympy.Symbol('kappa')
    alpha = sympy.Symbol('alpha')
    beta = sympy.Symbol('beta')
    m_3_2 = sympy.Symbol('m_3/2')

    # Part 1: Determine alpha^2 from the bosonic potential and Einstein's equations
    # The relation derived is R = 3 * kappa^2 * alpha^2
    alpha_squared_eq = sympy.Eq(R, 3 * kappa**2 * alpha**2)
    
    # Solve for alpha^2
    alpha_squared_sol = sympy.solve(alpha_squared_eq, alpha**2)[0]
    
    print("Step 1: Determine alpha^2")
    print("The relation between the Ricci scalar R and alpha is:")
    print(f"R = 3*kappa^2*alpha^2")
    print("Solving for alpha^2 gives:")
    # The final equation for alpha^2
    print("alpha**2 =")
    sympy.pprint(alpha_squared_sol)
    print("-" * 20)

    # Part 2: Determine beta from the gravitino mass terms
    # Relation 1: Mass from transformation law: m_3/2 = kappa*alpha/4
    # We assume alpha > 0.
    mass_from_transform_eq = sympy.Eq(m_3_2, kappa * alpha / 4)

    # Relation 2: Mass from Lagrangian comparison: alpha*kappa*beta = -m_3/2
    mass_from_lagrangian_eq = sympy.Eq(alpha * kappa * beta, -m_3_2)
    
    # Substitute m_3/2 from the first relation into the second
    final_beta_eq = mass_from_lagrangian_eq.subs(m_3_2, mass_from_transform_eq.rhs)
    
    # Solve for beta
    beta_sol = sympy.solve(final_beta_eq, beta)[0]

    print("Step 2: Determine beta")
    print("By comparing the gravitino mass terms from two different methods, we get the equation:")
    print(f"alpha*kappa*beta = -(kappa*alpha/4)")
    print("Solving for beta (assuming alpha and kappa are non-zero) gives:")
    # The final equation for beta
    print("beta =")
    sympy.pprint(beta_sol)
    print("-" * 20)

solve_sugra_parameters()