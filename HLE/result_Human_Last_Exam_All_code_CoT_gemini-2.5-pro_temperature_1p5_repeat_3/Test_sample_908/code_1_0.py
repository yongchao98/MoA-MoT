import sympy

def solve_electromagnetism_problem():
    """
    Solves for the total volume and surface charges in a cylindrical capacitor
    with non-uniform conductivity.
    """
    # Define symbolic variables
    V, eps, L, a, b, r, sigma0 = sympy.symbols('V epsilon L a b r sigma_0', real=True, positive=True)
    pi = sympy.pi
    I = sympy.Symbol('I', real=True, positive=True) # Total current

    # Step 1 & 2: Derive Electric Field E(r)
    # J_r = I / (2*pi*r*L)
    # sigma_r = sigma0 * r**2 / a**2
    # E_r = J_r / sigma_r
    E_r = (I / (2 * pi * r * L)) / (sigma0 * r**2 / a**2)

    # Step 3: Relate I to V
    V_expr = sympy.integrate(E_r, (r, a, b))
    # Solve for the current I in terms of the voltage V
    current_sol = sympy.solve(sympy.Eq(V, V_expr), I)
    I_val = current_sol[0]
    
    # Substitute I back into E_r to get the final expression for the electric field
    E_r_final = E_r.subs(I, I_val)
    E_r_final = sympy.simplify(E_r_final)

    # Step 4: Calculate volume charge density rho_v(r)
    # rho_v = eps * div(E) = eps * (1/r) * d(r*E_r)/dr
    div_E = (1/r) * sympy.diff(r * E_r_final, r)
    rho_v_r = eps * div_E

    # Step 5: Calculate total volume charge q_v
    # Integrate rho_v over the volume V' (2*pi*r*L*dr)
    q_v = sympy.integrate(rho_v_r * 2 * pi * r * L, (r, a, b))
    q_v = sympy.simplify(q_v)

    # Step 6: Calculate total surface charges q_s
    # Surface charge on inner electrode (r=a)
    # rho_s(a) = D_r(a) = eps * E_r(a)
    # q_s(a) = Area * rho_s(a) = (2*pi*a*L) * eps * E_r(a)
    q_s_a = sympy.simplify((2 * pi * a * L) * (eps * E_r_final.subs(r, a)))

    # Surface charge on outer electrode (r=b)
    # rho_s(b) = -D_r(b) = -eps * E_r(b)
    # q_s(b) = Area * rho_s(b) = (2*pi*b*L) * (-eps * E_r_final.subs(r, b))
    q_s_b = sympy.simplify((2 * pi * b * L) * (-eps * E_r_final.subs(r, b)))

    # Print the derived expressions for comparison with the answer choices.
    # The derived answers may look slightly different from the options due to
    # symbolic simplification, but are mathematically equivalent.
    print("Derived results:")
    print(f"Total volume charge = q_v = {q_v}")
    print(f"Total surface charge on inner electrode = q_s(r = a) = {q_s_a}")
    print(f"Total surface charge on outer electrode = q_s(r = b) = {q_s_b}")

    # Analysis of results:
    # My derived q_v is -4*pi*L*eps*V.
    # The derived surface charges are:
    # q_s(a) = 4*pi*L*eps*V / (1 - a**2/b**2)
    # q_s(b) = -4*pi*L*eps*V*a**2 / (b**2 * (1 - a**2/b**2))
    #
    # Comparing with the options, option D has the correct expressions for both q_s(a) and q_s(b).
    # However, the expression for q_v in option D is different and violates charge conservation.
    # Given that two out of three complex expressions match perfectly, option D is the most plausible
    # answer, likely containing a typo in the q_v term.

solve_electromagnetism_problem()