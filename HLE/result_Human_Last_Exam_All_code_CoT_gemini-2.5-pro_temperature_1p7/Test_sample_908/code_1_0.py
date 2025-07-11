import sympy as sp

def solve_electromagnetism_problem():
    """
    Solves for the total free volume and surface charges in a concentric
    cylindrical capacitor with inhomogeneous Ohmic material.
    """
    # 1. Define symbolic variables
    V, eps, L, a, b, r, sigma0 = sp.symbols('V epsilon L a b r sigma_0', real=True, positive=True)
    pi = sp.pi
    
    # 2. Define conductivity and find J_r
    # sigma(r) = sigma0 * r^2 / a^2
    # From continuity equation (div(J)=0), J_r = C/r
    C = sp.symbols('C') # Integration constant for current
    J_r = C / r

    # 3. Find E_r using Ohm's law
    sigma_r = sigma0 * r**2 / a**2
    E_r_expr = J_r / sigma_r
    # E_r_expr = (C/r) / (sigma0*r^2/a^2) = C * a^2 / (sigma0 * r^3)
    
    # 4. Find the constant C using the potential V
    # V = integral(E_r dr) from a to b
    V_integral = sp.integrate(E_r_expr, (r, a, b))
    # V_integral = C*a**2/sigma0 * integrate(r**-3, (r, a, b))
    # V_integral = C*a**2/sigma0 * [-1/(2*r**2)] from a to b
    # V_integral = C*a**2/(2*sigma0) * (1/a**2 - 1/b**2)
    # V_integral = C/(2*sigma0) * (1 - a**2/b**2)
    
    # Solve for C
    C_solved = sp.solve(sp.Eq(V, V_integral), C)[0]
    # C_solved = 2 * V * sigma0 / (1 - a**2/b**2)
    
    # Get the final expression for E_r
    E_r = E_r_expr.subs(C, C_solved)
    E_r = sp.simplify(E_r)
    # E_r = 2*V*a**2 / (r**3 * (1 - a**2/b**2))
    
    # 5. Calculate Volume Charge q_v
    # rho_v = eps * div(E) = eps * (1/r) * d/dr(r * E_r)
    r_E_r = r * E_r
    div_E = (1/r) * sp.diff(r_E_r, r)
    rho_v = eps * div_E
    
    # q_v = integral(rho_v * dVolume) = integral(rho_v * 2*pi*r*L*dr) from a to b
    integrand_qv = rho_v * 2 * pi * r * L
    q_v = sp.integrate(integrand_qv, (r, a, b))
    q_v = sp.simplify(q_v)

    # 6. Calculate Surface Charges q_s
    # rho_s(a) = eps * E_r(r=a)
    # q_s(a) = rho_s(a) * Area = eps * E_r(a) * 2*pi*a*L
    E_r_at_a = E_r.subs(r, a)
    q_s_a = eps * E_r_at_a * 2 * pi * a * L
    q_s_a = sp.simplify(q_s_a)
    
    # rho_s(b) = -eps * E_r(r=b)
    # q_s(b) = rho_s(b) * Area = -eps * E_r(b) * 2*pi*b*L
    E_r_at_b = E_r.subs(r, b)
    q_s_b = -eps * E_r_at_b * 2 * pi * b * L
    q_s_b = sp.simplify(q_s_b)

    # Print the results in a formatted way
    print("Derived Physical Quantities:")
    print("-" * 30)
    print("Total volume charge (q_v):")
    sp.pprint(q_v)
    print("\nTotal surface charge on inner electrode (q_s(a)):")
    sp.pprint(q_s_a)
    print("\nTotal surface charge on outer electrode (q_s(b)):")
    sp.pprint(q_s_b)

    # Final check: q_v + q_s_a + q_s_b should be 0
    total_charge = sp.simplify(q_v + q_s_a + q_s_b)
    print("\nVerification (Total Charge q_v + q_s(a) + q_s(b)):")
    sp.pprint(total_charge)


solve_electromagnetism_problem()
