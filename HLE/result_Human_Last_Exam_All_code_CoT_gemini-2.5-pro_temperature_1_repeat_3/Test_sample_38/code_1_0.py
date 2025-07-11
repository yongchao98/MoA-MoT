import sympy

def solve_gravity_mass():
    """
    Calculates the squared mass of the 6th degree of freedom in a modified
    linearized gravity theory.
    """

    # Step 1: Define the general framework.
    # A generic mass term for a spin-2 field h_munu can be written as:
    # L_mass = -(M^2 / 4) * (alpha * h_munu * h^munu + beta * h^2)
    # where h is the trace of h_munu.
    # This theory propagates a spin-2 particle (5 d.o.f.) and a spin-0 particle (1 d.o.f.).
    
    print("Step 1: State the general formulas for a generic graviton mass term.")
    
    # Define symbols for the general formula
    m_sq = sympy.Symbol('m^2')
    M_sq = sympy.Symbol('M^2')
    alpha = sympy.Symbol('alpha')
    beta = sympy.Symbol('beta')

    # According to the analysis by Boulware and Deser, the squared masses are:
    m2_sq_formula = M_sq * alpha / 2
    m0_sq_formula = M_sq * (alpha + 4*beta) / (2 * (alpha + 2*beta))

    print("The squared mass of the spin-2 part is given by: m_2^2 = M^2 * alpha / 2")
    print("The squared mass of the spin-0 part is given by: m_0^2 = M^2 * (alpha + 4*beta) / (2 * (alpha + 2*beta))\n")

    # Step 2: Identify parameters for the specific problem.
    # The given mass term is L_mass = -(m^2 / 2) * h_munu * h^munu.
    print("Step 2: Determine parameters for the Lagrangian L_mass = -(m^2 / 2) * h_munu * h^munu.")
    
    # Comparing this to the general form, the term proportional to h^2 is absent.
    beta_val = 0
    print(f"The h^2 term is absent, which implies: beta = {beta_val}")

    # Now we match the h_munu * h^munu part:
    # -(M^2 / 4) * alpha = -(m^2 / 2)
    # There is a freedom in choosing M^2 and alpha. A standard convention is to set alpha = 2.
    alpha_val = 2
    print(f"We adopt the standard convention: alpha = {alpha_val}")
    
    # With alpha = 2, the equation for the mass scales becomes:
    # -(M^2 / 4) * 2 = -(m^2 / 2)  =>  -M^2 / 2 = -m^2 / 2
    M_sq_val_sym = m_sq
    print(f"Solving for M^2 gives: M^2 = m^2\n")

    # Step 3: Calculate the squared mass of the sixth degree of freedom (the spin-0 part).
    print("Step 3: Calculate the squared masses of the propagating particles.")
    
    # As a sanity check, we first verify the spin-2 mass given in the problem.
    m2_sq_final = m2_sq_formula.subs({M_sq: M_sq_val_sym, alpha: alpha_val})
    print(f"Sanity Check: Calculating the spin-2 mass-squared:")
    print(f"m_2^2 = M^2 * alpha / 2 = m^2 * {alpha_val} / 2 = {m2_sq_final}")
    print("This confirms the problem statement that 5 d.o.f. have mass-squared m^2.\n")

    # Now, we perform the main calculation for the spin-0 mass.
    print("Calculating the squared mass of the sixth (spin-0) degree of freedom:")
    
    # Substitute the determined parameter values into the general formula for m_0^2
    num_val = alpha_val + 4 * beta_val
    den_val = 2 * (alpha_val + 2 * beta_val)
    m0_sq_final = m0_sq_formula.subs({M_sq: M_sq_val_sym, alpha: alpha_val, beta: beta_val})
    
    print(f"Using the formula: m_0^2 = M^2 * (alpha + 4*beta) / (2 * (alpha + 2*beta))")
    print(f"Substituting M^2=m^2, alpha={alpha_val}, beta={beta_val}:")
    print(f"m_0^2 = m^2 * ({alpha_val} + 4*{beta_val}) / (2 * ({alpha_val} + 2*{beta_val}))")
    print(f"m_0^2 = m^2 * ({num_val}) / ({den_val})")
    print(f"The result is: m_0^2 = {m0_sq_final}")

solve_gravity_mass()