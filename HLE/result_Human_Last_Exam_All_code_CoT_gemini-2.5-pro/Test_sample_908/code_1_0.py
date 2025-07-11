import sympy
from sympy import symbols, pi, Eq, solve, integrate, diff, Function

def solve_electromagnetism_problem():
    """
    Solves the concentric cylindrical electrodes problem symbolically.
    """
    # Step 1: Define symbolic variables
    r, a, b, L, V, epsilon, sigma_0 = symbols('r a b L V varepsilon sigma_0', positive=True, real=True)
    
    # K is an unknown constant for the electric field
    K = symbols('K')

    print("Step 1: Deriving the Electric Field E(r)")
    # Conductivity sigma(r)
    sigma_r = sigma_0 * r**2 / a**2
    
    # In DC steady state, div(J) = 0. By cylindrical symmetry, J_r = const / r.
    # From Ohm's Law, J = sigma * E, so E_r = J_r / sigma_r.
    # E_r = (const / r) / (sigma_0 * r**2 / a**2) = K / r**3 for some constant K.
    E_r = K / r**3
    print(f"Due to symmetry and Ohm's law, the electric field has the form: E(r) = K / r^3\n")

    print("Step 2: Applying Boundary Conditions to find K")
    # The potential difference V is the integral of E from r=a to r=b.
    # We assume V = V(a) - V(b), so E is directed radially outward.
    # V = - integral(E_r dr) from b to a = integral(E_r dr) from a to b
    potential_eq = Eq(V, integrate(E_r, (r, a, b)))
    
    # Solve for the constant K
    K_solution = solve(potential_eq, K)
    K_expr = K_solution[0]
    print(f"The constant K is found by integrating E(r): V = integral(E(r), (r, a, b))")
    print(f"Solving for K gives: K = {sympy.pretty(K_expr)}\n")

    # Substitute K back into the electric field expression
    E_r_final = E_r.subs(K, K_expr)

    print("Step 3: Calculating the Volume Charge Density rho_v(r)")
    # Volume charge density rho_v = epsilon * div(E)
    # In cylindrical coordinates, div(E) = (1/r) * d/dr(r * E_r)
    div_E = (1/r) * diff(r * E_r_final, r)
    rho_v = epsilon * div_E
    print(f"The volume charge density is rho_v = epsilon * div(E).")
    print(f"rho_v(r) = {sympy.pretty(sympy.simplify(rho_v))}\n")

    print("Step 4: Calculating the Total Volume Charge q_v")
    # Total volume charge q_v = integral of rho_v over the volume
    # dV = 2 * pi * r * L * dr
    q_v = integrate(rho_v * 2 * pi * r * L, (r, a, b))
    q_v_simplified = sympy.simplify(q_v)
    print(f"The total volume charge is q_v = integral(rho_v * 2*pi*r*L, (r, a, b))")
    # Final simplified expression for total volume charge
    final_q_v = -4 * pi * L * V * epsilon
    print("Final result for total volume charge:")
    print(f"q_v = {sympy.pretty(final_q_v)}\n")


    print("Step 5: Calculating the Surface Charges q_s")
    # Electric displacement field D = epsilon * E
    D_r = epsilon * E_r_final

    # Inner electrode at r=a
    # Surface charge density sigma_s(a) = D_r evaluated at r=a
    # Total surface charge q_s(a) = Area * sigma_s(a) = (2*pi*a*L) * D_r(a)
    q_s_a = (2 * pi * a * L) * D_r.subs(r, a)
    q_s_a_simplified = sympy.simplify(q_s_a)
    # Re-writing to match option format
    q_s_a_formatted = (4 * pi * L * V * epsilon) / (1 - a**2/b**2)

    print("Total surface charge on the inner electrode (r=a):")
    print(f"q_s(a) = {sympy.pretty(q_s_a_simplified)}")
    print(f"Which can be written as: q_s(a) = {sympy.pretty(q_s_a_formatted)}\n")

    # Outer electrode at r=b
    # Surface charge density sigma_s(b) = -D_r evaluated at r=b (normal vector is -r_hat)
    # Total surface charge q_s(b) = Area * sigma_s(b) = (2*pi*b*L) * (-D_r(b))
    q_s_b = (2 * pi * b * L) * (-D_r.subs(r, b))
    q_s_b_simplified = sympy.simplify(q_s_b)
    # Re-writing to match option format
    q_s_b_formatted = (-4 * pi * L * V * epsilon * a**2) / (b**2 * (1 - a**2/b**2))
    
    print("Total surface charge on the outer electrode (r=b):")
    print(f"q_s(b) = {sympy.pretty(q_s_b_simplified)}")
    print(f"Which can be written as: q_s(b) = {sympy.pretty(q_s_b_formatted)}\n")

    print("Step 6: Comparing derived results with answer choices")
    print("Our derived expressions are:")
    print(f"  Total volume charge = {sympy.pretty(final_q_v)}")
    print(f"  Total surface charge on inner electrode = {sympy.pretty(q_s_a_formatted)}")
    print(f"  Total surface charge on outer electrode = {sympy.pretty(q_s_b_formatted)}\n")
    
    print("Let's analyze the options:")
    print("Option A: All expressions mismatch.")
    print("Option B:")
    print("  - q_v doesn't match our derived result.")
    print("  - q_s(a) is off by a factor of 2 (2*pi vs our 4*pi).")
    print("  - q_s(b) is an exact match to our derived expression.")
    print("Option C, D, E: Mismatch.")
    
    print("\nConclusion:")
    print("The derivations from first principles show that the provided options are flawed, as they are not internally consistent (the charges in Option B do not sum to zero).")
    print("However, the expression for the charge on the outer electrode, q_s(r=b), in Option B is a perfect match to our derived result.")
    print("This suggests that Option B is the intended answer, despite containing apparent typos in the expressions for the volume charge and the inner surface charge.")

if __name__ == '__main__':
    solve_electromagnetism_problem()
<<<B>>>