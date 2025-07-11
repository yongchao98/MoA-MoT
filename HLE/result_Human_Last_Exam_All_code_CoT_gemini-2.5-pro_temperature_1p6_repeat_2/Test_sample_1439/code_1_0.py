def solve_critical_exponent_order():
    """
    This script explains and demonstrates the order in the coupling constant 'u'
    at which the critical exponent nu acquires its first non-vanishing contribution.
    """

    # We use string representations for the symbolic variables to make the printout clear.
    u = 'u'
    N_components = 1  # For the N=1 universality class (like the Ising model).

    # Step 1: State the fundamental relationship from Renormalization Group theory.
    print("Within the Renormalization Group (RG) framework, the critical exponent nu is related to")
    print("the anomalous dimension gamma_r at the non-trivial fixed point.")
    print("The relationship is expressed as a series in the coupling constant u:")
    print("  1/nu = 2 + gamma_r(u)")
    print("-" * 50)

    # Step 2: Provide the one-loop (first order in u) expansion for gamma_r(u).
    # For a theory with N components, gamma_r(u) = -((N+2)/6) * u + O(u^2)
    gamma_r_numerator = N_components + 2
    gamma_r_denominator = 6

    print(f"The one-loop calculation for a theory with N={N_components} component(s) gives:")
    print(f"  gamma_r(u) = -({gamma_r_numerator}/{gamma_r_denominator}) * {u}^1 + O({u}^2)")
    # For N=1, this simplifies to -3/6 * u = -1/2 * u
    simplified_gamma_r_num = 1
    simplified_gamma_r_den = 2
    print(f"  which simplifies to: -({simplified_gamma_r_num}/{simplified_gamma_r_den}) * {u}^1")
    print("-" * 50)


    # Step 3: Substitute this into the equation for 1/nu.
    print("Substituting this into the equation for nu gives the expansion for 1/nu:")
    print(f"  1/nu = 2 - ({simplified_gamma_r_num}/{simplified_gamma_r_den}) * {u}^1 + O({u}^2)")
    print("-" * 50)

    # Step 4: Invert the series to find the expansion for nu itself.
    # nu ≈ (1/2) * [1 - (1/4)*u]^-1 ≈ (1/2) * (1 + (1/4)*u) = 1/2 + 1/8 * u
    mean_field_num = 1
    mean_field_den = 2
    correction_num = 1
    correction_den = 8
    order_of_u = 1

    print("To find nu, we invert the expression. The mean-field value (at u=0) is 1/2.")
    print("The full expansion for nu starts as:")
    print(f"  nu(u) = ({mean_field_num}/{mean_field_den}) + ({correction_num}/{correction_den}) * {u}^{order_of_u} + O({u}^2)")
    print("-" * 50)


    # Step 5: Conclude the order of the first contribution.
    print("The first non-vanishing contribution to nu beyond its mean-field value (1/2)")
    print(f"is the term '({correction_num}/{correction_den}) * u^{order_of_u}'.")
    print(f"\nThis contribution appears at order {order_of_u} in the coupling constant u.")

solve_critical_exponent_order()