import sympy
from sympy import symbols, pi, integrate, diff, Eq, solve, pretty_print

def solve_electromagnetism_problem():
    """
    This function symbolically derives the charges for the given concentric
    cylindrical electrode problem using the sympy library.
    """
    # Define symbolic variables
    r, a, b, L, V, epsilon, sigma0 = symbols('r a b L V epsilon sigma_0', positive=True, real=True)
    K = symbols('K', real=True)  # Proportionality constant for current density

    # --- Step 1: Define E-field in terms of an unknown constant K ---
    # From ∇ ⋅ J = 0, we get J_r = K/r
    J_r = K / r
    # Given conductivity σ(r)
    sigma = sigma0 * r**2 / a**2
    # From Ohm's Law, E_r = J_r / σ
    E_r = J_r / sigma

    # --- Step 2: Use Voltage V to solve for K ---
    # V = ∫ E_r dr from a to b
    voltage_eq = Eq(V, integrate(E_r, (r, a, b)))
    # Solve for the constant K
    K_solution = solve(voltage_eq, K)
    if not K_solution:
        print("Could not solve for K.")
        return
    K_val = K_solution[0]

    # Substitute K back into the expression for the electric field
    E_r_final = E_r.subs(K, K_val)

    # --- Step 3: Calculate Volume Charge q_v ---
    # ρ_v = ε * ∇ ⋅ E
    div_E = (1/r) * diff(r * E_r_final, r)
    rho_v = epsilon * div_E
    # q_v = ∫ ρ_v d(Volume)
    q_v = integrate(rho_v * 2 * pi * r * L, (r, a, b))

    # --- Step 4: Calculate Surface Charges q_s(a) and q_s(b) ---
    # q_s(a) = Area * ε * E_r(a)
    E_r_at_a = E_r_final.subs(r, a)
    q_s_a = (2 * pi * a * L * epsilon * E_r_at_a)

    # q_s(b) = Area * (-ε) * E_r(b)
    E_r_at_b = E_r_final.subs(r, b)
    q_s_b = (2 * pi * b * L * (-epsilon) * E_r_at_b)

    # --- Step 5: Print the derived expressions ---
    print("The derived expressions for the charges are:")
    print("\nTotal volume charge (q_v):")
    pretty_print(sympy.simplify(q_v))

    print("\nTotal surface charge on inner electrode q_s(r=a):")
    pretty_print(sympy.simplify(q_s_a))

    print("\nTotal surface charge on outer electrode q_s(r=b):")
    pretty_print(sympy.simplify(q_s_b))

    # Compare with Option B
    q_s_b_option_B = -4 * pi * L * V * epsilon * a**2 / (b**2 * (1 - a**2/b**2))
    print("\nComparing derived q_s(b) with Option B:")
    print("Is derived q_s(b) equal to Option B's q_s(b)?", sympy.simplify(q_s_b - q_s_b_option_B) == 0)


solve_electromagnetism_problem()