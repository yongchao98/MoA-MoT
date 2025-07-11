def solve_critical_exponent_order():
    """
    This function determines and prints the order of the first non-vanishing contribution
    to the critical exponent ν from the coupling constant u in φ⁴ theory.
    """

    # The final expanded equation for the critical exponent ν as a function of the coupling u
    # has the form: ν(u) = ν₀ + C * u^k + O(u²)
    # We identify the numbers in this equation based on the perturbative analysis.

    # 1. The mean-field value (ν at u=0)
    nu_mean_field = 0.5

    # 2. The power 'k' of the coupling 'u' in the first non-zero correction term.
    # This comes from the one-loop calculation.
    power_of_u_in_first_correction = 1

    # 3. The order of the next term in the expansion.
    order_of_next_term = 2

    # A string to represent the constant coefficient of the first correction term.
    constant_C_string = "(c₁/4)"

    print("The critical exponent ν can be expressed as a power series in the coupling constant u:")
    print(f"ν(u) = {nu_mean_field} + {constant_C_string} * u^{power_of_u_in_first_correction} + O(u^{order_of_next_term})")
    print("-" * 50)
    print("This equation shows that the initial non-vanishing contribution to ν (the term modifying the mean-field value of 0.5) is proportional to u raised to the first power.")
    print("\nTherefore, the order of this contribution is:")
    print(power_of_u_in_first_correction)

solve_critical_exponent_order()