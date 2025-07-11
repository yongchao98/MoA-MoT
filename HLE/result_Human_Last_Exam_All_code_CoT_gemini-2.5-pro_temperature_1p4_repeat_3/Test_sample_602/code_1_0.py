import math

def calculate_l(n: int):
    """
    Calculates the exact value of the function l(n) for n >= 5.

    The formula is derived step-by-step as outlined in the plan.
    The final simplified expression for l(n) is:
    l(n) = 2 + 2/n^2 - (4*n - 2)/n^2 * sqrt(n^2 - n + 1)
    """
    if n < 5:
        raise ValueError("The function is defined for n >= 5.")

    # Calculate the components of the formula
    term1 = 2
    term2_numerator = 2
    term2_denominator = n**2
    
    term3_sqrt_arg = n**2 - n + 1
    term3_sqrt_val = math.sqrt(term3_sqrt_arg)
    
    term3_numerator_coeff = -(4 * n - 2)
    term3_denominator = n**2
    
    # Calculate the final value
    result = term1 + (term2_numerator / term2_denominator) + (term3_numerator_coeff / term3_denominator) * term3_sqrt_val

    # Print the equation with its numerical components
    print(f"For n = {n}, the equation for l(n) is:")
    print(f"l(n) = c1 + c2/n^2 + c3/n^2 * sqrt(n^2 - n + 1)")
    print(f"where:")
    print(f"c1 = {term1}")
    print(f"c2 = {term2_numerator}")
    print(f"c3 = {term3_numerator_coeff}")
    print(f"n^2 - n + 1 = {term3_sqrt_arg}")
    
    print("\nThe final calculated value is:")
    print(f"l({n}) = {result}")

# Example calculation for n=5
calculate_l(5)