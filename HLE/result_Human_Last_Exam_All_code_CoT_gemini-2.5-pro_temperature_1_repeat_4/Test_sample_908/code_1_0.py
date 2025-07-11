import sympy
from sympy import symbols, pi, Function, dsolve, Eq, solve, integrate, Derivative

def solve_electrodynamics_problem():
    """
    Solves the concentric cylinder electrodynamics problem symbolically.
    """
    # Define symbols
    r, a, b, L, V, epsilon, sigma_0 = symbols('r a b L V varepsilon sigma_0', real=True, positive=True)
    
    # --- Step 1 & 2: Find Electric Field E(r) ---
    # From ∇ ⋅ J = 0 in cylindrical coordinates, r*J_r is constant. J = σE.
    # So r * σ(r) * E(r) = constant.
    # Let's find the total current I first.
    # E_r(r) = I / (2*pi*r*L * σ(r))
    sigma = sigma_0 * r**2 / a**2
    
    # We need to find I in terms of V.
    # V = integral from a to b of E_r(r) dr
    I = symbols('I') # Total Current
    E_r_expr = I / (2 * pi * r * L * sigma)
    
    # V = ∫[a,b] E_r dr
    voltage_eq = Eq(V, integrate(E_r_expr, (r, a, b)))
    
    # Solve for I
    # integrate(1/r^3) = -1/(2*r^2)
    # integral from a to b is -1/(2*b^2) - (-1/(2*a^2)) = 1/2 * (1/a^2 - 1/b^2)
    # V = I*a**2 / (2*pi*L*sigma_0) * (1/2 * (1/a**2 - 1/b**2))
    # V = I / (4*pi*L*sigma_0) * (1 - a**2/b**2)
    I_solved = solve(voltage_eq, I)[0]
    
    # Substitute I back into E_r expression
    E_r = E_r_expr.subs(I, I_solved).simplify()

    # --- Step 3 & 4: Find Total Volume Charge q_v ---
    # ρ_v = ε * ∇⋅E
    # ∇⋅E = (1/r) * d(r*E_r)/dr
    div_E = (1/r) * Derivative(r * E_r, r).doit().simplify()
    rho_v = epsilon * div_E
    
    # q_v = ∫ ρ_v d(Volume) = ∫[0,L] ∫[0,2π] ∫[a,b] ρ_v * r dr dφ dz
    q_v = integrate(rho_v * r * 2 * pi * L, (r, a, b)).simplify()

    # --- Step 5 & 6: Find Total Surface Charges q_s ---
    # Inner electrode at r=a
    # σ_s(a) = D_r(a) = ε * E_r(a)
    # Area = 2*pi*a*L
    q_s_a = (2 * pi * a * L * epsilon * E_r.subs(r, a)).simplify()
    
    # Outer electrode at r=b
    # The normal vector points into the conductor, so σ_s(b) = -D_r(b) = -ε * E_r(b)
    # Area = 2*pi*b*L
    q_s_b = (-1 * 2 * pi * b * L * epsilon * E_r.subs(r, b)).simplify()

    # --- Step 7: Compare with options ---
    # The derived results are:
    # q_v = -4*pi*L*V*epsilon
    # q_s(a) = 4*pi*L*V*epsilon / (1 - a**2/b**2)
    # q_s(b) = -4*pi*L*V*epsilon*a**2 / (b**2 * (1 - a**2/b**2))
    #
    # A self-consistency check is q_v + q_s(a) + q_s(b) = 0.
    # Let's check: (-4) + (4 / (1-x)) + (-4x / (1-x)) = (-4(1-x) + 4 - 4x) / (1-x)
    # = (-4 + 4x + 4 - 4x) / (1-x) = 0. The results are consistent.
    #
    # None of the options perfectly match these results. The options provided seem to be flawed.
    # Let's analyze Option B, as it's the most plausible choice despite its errors.
    
    denom = (1 - a**2/b**2)
    
    # Expressions from Option B
    q_v_B = -4 * V * epsilon * pi * L / denom
    q_s_a_B = 2 * pi * L * V * epsilon / denom
    q_s_b_B = -4 * pi * L * V * epsilon * a**2 / (b**2 * denom)

    print("--- Derived Results ---")
    print(f"Total volume charge (q_v) = {q_v}")
    print(f"Surface charge on inner electrode (q_s(a)) = {q_s_a}")
    print(f"Surface charge on outer electrode (q_s(b)) = {q_s_b}")
    print("\n--- Expressions from Option B ---")
    print(f"Total volume charge (q_v) = {q_v_B}")
    print(f"Surface charge on inner electrode (q_s(a)) = {q_s_a_B}")
    print(f"Surface charge on outer electrode (q_s(b)) = {q_s_b_B}")
    
    print("\n--- Analysis ---")
    print("The derived results are self-consistent (sum to zero), but do not perfectly match any option.")
    print("Comparing with Option B:")
    print(f" - q_s(b) matches the derived result: {q_s_b == q_s_b_B}")
    print(f" - q_s(a) differs from the derived result by a factor of 2.")
    print(f" - q_v has the correct sign but an incorrect functional form compared to the derived result.")
    print("Given the inconsistencies in all options, Option B is the most likely intended answer, assuming typos in its expressions for q_v and q_s(a).")


solve_electrodynamics_problem()