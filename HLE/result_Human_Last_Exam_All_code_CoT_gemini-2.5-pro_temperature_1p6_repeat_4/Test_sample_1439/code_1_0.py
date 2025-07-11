import math

def print_critical_exponent_expansion():
    """
    This function explains and prints the perturbative expansion for the critical
    exponent ν in φ⁴ theory.
    """

    # Introduction to the context
    print("In the renormalization group analysis of φ⁴ theory, the critical exponent ν is expanded as a power series in the coupling constant u.")
    print("The expansion around the mean-field value is:")
    print("ν(u) = ν₀ + c₁*u¹ + c₂*u² + ...")
    print("-" * 30)

    # Define the terms in the equation
    nu_0_numerator = 1
    nu_0_denominator = 2
    
    first_order_coeff_symbol = "c₁"
    first_order_variable = "u"
    first_order_power = 1

    higher_order_symbol = "O(u²)"

    # Print the explanation of each term
    print(f"The zeroth-order term (mean-field value) is ν₀:")
    print(f"ν₀ = {nu_0_numerator}/{nu_0_denominator}")
    print("\nThis value corresponds to the non-interacting theory (u=0).")
    print("-" * 30)

    print(f"The first correction term at the non-trivial fixed point is:")
    print(f"Contribution = {first_order_coeff_symbol} * {first_order_variable}^{first_order_power}")
    print("\nOne-loop calculations show that the coefficient c₁ is non-zero.")
    print("This is the initial non-vanishing contribution to ν beyond its mean-field value.")
    print("-" * 30)

    # Print the final, full equation
    print("The full equation showing each number and symbol for the expansion is:")
    print(f"ν = {nu_0_numerator}/{nu_0_denominator} + {first_order_coeff_symbol} * {first_order_variable}^{first_order_power} + {higher_order_symbol}")
    print("\nConclusion: The critical exponent ν acquires its initial non-vanishing contribution at the first order (power of 1) in the coupling constant u.")

if __name__ == "__main__":
    print_critical_exponent_expansion()