import sympy
from sympy import symbols, Eq, solve, integrate, diff, pi

def solve_electromagnetics_problem():
    """
    Solves for the volume and surface charges in a concentric cylindrical capacitor
    with non-uniform conductivity using symbolic mathematics.
    """
    # Define the symbolic variables
    V, epsilon, L, a, b, r, sigma0 = symbols('V epsilon L a b r sigma_0', real=True, positive=True)
    I = symbols('I', real=True) # Total DC current

    # Step 1 & 2: Find J(r) and E(r)
    # From ∇⋅J=0 in cylindrical coords -> J(r) = I / (2*pi*r*L)
    # Given sigma(r) = sigma0 * r^2 / a^2
    # From Ohm's Law, E = J/sigma
    J_r = I / (2 * pi * r * L)
    sigma_r = sigma0 * r**2 / a**2
    E_r = J_r / sigma_r

    # Step 3: Relate I to V
    # V = integral(E_r dr) from a to b
    V_expr = integrate(E_r, (r, a, b))
    
    # Solve for I in terms of V
    I_solution = solve(Eq(V, V_expr), I)[0]

    # Substitute I back into E_r to get E_r(V)
    E_r_final = E_r.subs(I, I_solution)

    # Step 4: Find the volume charge density rho_f
    # rho_f = epsilon * ∇⋅E. For radial E, ∇⋅E = (1/r) * d(r*E_r)/dr
    div_E = (1/r) * diff(r * E_r_final, r)
    rho_f = epsilon * div_E

    # Step 5: Calculate total volume charge q_v
    # q_v = integral(rho_f * dVolume) where dVolume = 2*pi*r*L*dr
    q_v = integrate(rho_f * (2 * pi * r * L), (r, a, b))

    # Step 6: Calculate total surface charges q_s(a) and q_s(b)
    # q_s(a) = Area(a) * D_r(a) = (2*pi*a*L) * epsilon * E_r(r=a)
    q_s_a = (2 * pi * a * L) * epsilon * E_r_final.subs(r, a)
    
    # q_s(b) = Area(b) * (-D_r(b)) = (2*pi*b*L) * (-epsilon) * E_r(r=b)
    q_s_b = (2 * pi * b * L) * (-epsilon) * E_r_final.subs(r, b)

    # Simplify expressions for printing
    q_v_final = sympy.simplify(q_v)
    q_s_a_final = sympy.simplify(q_s_a)
    q_s_b_final = sympy.simplify(q_s_b)

    # Re-format for comparison with options
    q_s_a_final_reformat = (4 * pi * L * V * epsilon) / (1 - (a**2/b**2))
    q_s_b_final_reformat = (-4 * pi * L * V * epsilon * a**2) / (b**2 * (1-a**2/b**2))

    print("Derived Total Volume Charge (q_v):")
    sympy.pprint(q_v_final)
    
    print("\nDerived Total Surface Charge on Inner Electrode q_s(r=a):")
    sympy.pprint(q_s_a_final_reformat)

    print("\nDerived Total Surface Charge on Outer Electrode q_s(r=b):")
    sympy.pprint(q_s_b_final_reformat)

    # Check for charge neutrality
    neutrality_check = sympy.simplify(q_v_final + q_s_a_final + q_s_b_final)
    print(f"\nNeutrality Check (q_v + q_s(a) + q_s(b)): {neutrality_check}")

solve_electromagnetics_problem()
