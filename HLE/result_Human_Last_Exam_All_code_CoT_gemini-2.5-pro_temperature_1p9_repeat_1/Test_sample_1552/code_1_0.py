import sympy

def find_heat_kernel_coefficient():
    """
    Calculates the density of the second coefficient (a_2) in the heat kernel expansion
    for a massless gauged Dirac spinor field in 4D spacetime.
    """

    # Define symbols
    R, N_c, pi = sympy.symbols('R N_c pi')
    I_spin, I_gauge = sympy.symbols('I_spin I_gauge')

    # Step 1: The Operator and Heat Kernel expansion
    # The operator is H = D^2, the square of the Dirac operator.
    # We want the coefficient a_1 in the expansion Tr(exp(-tH)) ~ (4*pi*t)^(-2) * sum(a_n * t^n)
    # This coefficient is often called a_2.

    # Step 2: Lichnerowicz Formula and relation to standard heat kernel form
    # The Lichnerowicz formula gives D^2 = Delta_B + E_H, where Delta_B is the Bochner Laplacian,
    # and E_H = (1/4)*R*I + (1/2)*sigma_munu*F_munu.
    # The standard formula for the a_1 coefficient is for an operator L = -Delta_B + E_L.
    # a_1(L) = Tr( (1/6)*R*I + E_L ).
    # Relating H to L, we have E_L = -E_H.
    # So, a_1(H) = Tr( (1/6)*R*I - E_H ).
    
    # Substituting E_H:
    # a_1_integrand = Tr( (1/6)*R*I - ( (1/4)*R*I + (1/2)*sigma_munu*F_munu ) )
    # a_1_integrand = Tr( (1/6 - 1/4)*R*I - (1/2)*sigma_munu*F_munu )
    
    coeff_R_num = 1
    coeff_R_den_1 = 6
    coeff_R_den_2 = 4

    term_R_factor = sympy.Rational(coeff_R_num, coeff_R_den_1) - sympy.Rational(coeff_R_num, coeff_R_den_2)

    # Step 3: Compute the trace
    # The trace of the second term involving sigma_munu (Lorentz generators) is zero
    # because Tr_spin(sigma_munu) = 0.
    # So we only need to compute the trace of the R-term.
    # Tr(c*R*I) = c * R * Tr(I_spin x I_gauge) = c * R * Tr_spin(I_spin) * Tr_gauge(I_gauge)

    dim_spinor = 4  # Dimension of Dirac spinors in 4D
    dim_gauge = N_c # Dimension of the gauge group representation

    traced_term = term_R_factor * R * dim_spinor * dim_gauge

    # Step 4: Apply the normalization factor
    # The full density is (1/((4*pi)^(d/2))) * a_1_integrand, with d=4.
    # Normalization factor is 1 / (16 * pi**2)
    norm_factor_den = 16

    a2_density = traced_term / (norm_factor_den * pi**2)
    
    # --- Output Results ---
    print("The density of the second heat kernel coefficient (a_2) is calculated as follows:")
    print("-" * 70)
    
    print("1. The general formula for the untraced integrand is:")
    print(f"   b_2(x) = (1/6 - 1/4)*R*I - (1/2)*sigma_munu*F_munu")
    print(f"   Which simplifies to: ({term_R_factor})*R*I, as Tr(sigma_munu*F_munu) = 0.\n")

    print("2. Taking the trace over spinor (dim=4) and gauge (dim=N_c) spaces:")
    calc_step_2 = f"   Tr[({term_R_factor})*R*I] = ({term_R_factor}) * R * {dim_spinor} * N_c = {sympy.simplify(term_R_factor * dim_spinor)} * N_c * R"
    print(calc_step_2 + "\n")

    print("3. Applying the normalization factor 1/(16*pi^2):")
    final_num = sympy.simplify(term_R_factor * dim_spinor).p
    final_den_part1 = sympy.simplify(term_R_factor * dim_spinor).q
    total_den = final_den_part1 * norm_factor_den
    print(f"   Density = ({final_num} * N_c * R) / ({final_den_part1} * {norm_factor_den} * pi^2)\n")

    print("Final expression for the coefficient density b_2(x):")
    final_expr = sympy.simplify(a2_density)
    
    # Explicitly print the numbers in the final simplified equation
    num_final, den_final = final_expr.as_numer_den()
    print("b_2(x) = ", end="")
    sympy.pprint(final_expr)
    
    print("\n--- Equation Breakdown ---")
    num_parts = num_final.as_ordered_factors()
    den_parts = den_final.as_ordered_factors()
    
    num_const = [n for n in num_parts if n.is_Number]
    den_const = [d for d in den_parts if d.is_Number]
    
    print(f"Numerical constant: {num_const[0] if num_const else 1}/{den_const[0] if den_const else 1}")
    print(f"Symbolic part: N_c * R / pi**2")


find_heat_kernel_coefficient()
<<<(-N_c*R)/(12*pi**2)>>>