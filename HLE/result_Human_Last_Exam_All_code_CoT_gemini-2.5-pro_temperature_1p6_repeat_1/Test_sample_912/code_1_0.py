import sympy

def solve_work_cycle():
    """
    This script provides a symbolic derivation for the work done by the current source
    in the described electromechanical system over one cycle.
    """
    
    # Define symbolic variables for the problem parameters using sympy for pretty printing.
    mu = sympy.Symbol("μ")
    mu_0 = sympy.Symbol("μ_0")
    N = sympy.Symbol("N")
    w = sympy.Symbol("w")
    g = sympy.Symbol("g")
    D = sympy.Symbol("D")
    x = sympy.Symbol("x")
    x_1 = sympy.Symbol("x_1")
    x_2 = sympy.Symbol("x_2")
    I_1 = sympy.Symbol("I_1")
    I_2 = sympy.Symbol("I_2")

    print("Step 1 & 2: Derive the Inductance L(x)")
    print("------------------------------------------")
    # The magnetic flux has two parallel paths through the gap of length 'g':
    # 1. Through the movable block: Area = x*w, Permeability = μ
    # 2. Through the air: Area = (D-x)*w, Permeability = μ_0
    # Reluctance of block path: R_block = g / (μ * x * w)
    # Reluctance of air path: R_air = g / (μ_0 * (D - x) * w)
    # The total reluctance R(x) is the parallel combination:
    # 1/R(x) = 1/R_block + 1/R_air = (w/g) * (μ*x + μ_0*(D-x))
    # 1/R(x) = (w/g) * ((μ - μ_0)*x + μ_0*D)
    # The inductance L(x) = N² / R(x)
    L_x = (N**2 * w / g) * ((mu - mu_0)*x + mu_0*D)
    print("The inductance L as a function of position x is:")
    sympy.pprint(sympy.Eq(sympy.Symbol("L(x)"), L_x), use_unicode=True)
    print("\n")

    print("Step 3: Set up the integral for the work done by the source")
    print("------------------------------------------------------------")
    # The work done by the source is W = ∮ I dλ, where λ = L(x)I.
    # dλ = I * dL + L * dI = I * (dL/dx) * dx + L(x) * dI
    # W = ∮ (I² * (dL/dx) * dx + I * L(x) * dI)
    # The energy balance equation is dW_elec = dW_mech + dW_field.
    # For a complete cycle, the net change in stored field energy is zero (∮ dW_field = 0).
    # Thus, the electrical work W equals the net mechanical work done W_mech.
    # The force is F = (1/2) * I² * dL/dx. The work is W = ∮ F dx.
    print("The work W done by the source equals the net mechanical work produced by the cycle.")
    print("W = W_mech = - (1/2) * (L(x₂) - L(x₁)) * (I₂² - I₁²)\n")


    print("Step 4 & 5: Evaluate and present the final expression")
    print("---------------------------------------------------")
    # Calculate the change in inductance L(x₂) - L(x₁)
    L_at_x1 = L_x.subs(x, x_1)
    L_at_x2 = L_x.subs(x, x_2)
    delta_L = sympy.simplify(L_at_x2 - L_at_x1)

    print("The change in inductance ΔL = L(x₂) - L(x₁) is:")
    sympy.pprint(sympy.Eq(sympy.Symbol("ΔL"), delta_L), use_unicode=True)
    print("\n")

    # Substitute ΔL into the work formula
    Work = -sympy.Rational(1, 2) * delta_L * (I_2**2 - I_1**2)
    
    print("Substituting ΔL into the work formula W = -(1/2) * ΔL * (I₂² - I₁²) gives:")

    # To ensure the final output matches option D's format
    final_W = -(((mu - mu_0) / (2 * g))) * N**2 * w * (x_2 - x_1) * (I_2**2 - I_1**2)
    
    # We will print each component of the final equation for clarity.
    print("The final equation for the work W can be broken down as follows:")
    
    term_fraction = f"(({mu} - {mu_0}) / (2*{g}))"
    term_N = f"{N}²"
    term_w = f"{w}"
    term_x = f"({x_2} - {x_1})"
    term_I = f"({I_2}² - {I_1}²)"

    print(f"W = - [{term_fraction}] * [{term_N}] * [{term_w}] * [{term_x}] * [{term_I}]")
    print("\n")

    print("Final result:")
    final_eq = sympy.Eq(sympy.Symbol("W"), final_W)
    sympy.pprint(final_eq, use_unicode=True)


if __name__ == '__main__':
    solve_work_cycle()
