import math

def calculate_l(d):
    """
    Calculates the value of l(d).
    
    The problem simplifies to finding the minimum of the Busemann function,
    which results in the expression: l(d) = ln((sqrt(d)-1)/(sqrt(d)+1)).
    
    Args:
        d (int): An integer d >= 2.
    
    Returns:
        float: The calculated value of l(d).
    """
    if d < 2:
        raise ValueError("d must be greater than or equal to 2.")
    
    sqrt_d = math.sqrt(d)
    numerator = sqrt_d - 1
    denominator = sqrt_d + 1
    
    # Check for non-positive argument to log
    if numerator <= 0 or denominator <= 0:
      # This case happens only if d<=1 which is checked above.
      # Return NaN or raise an error for robustness.
      return float('nan')
      
    value = math.log(numerator / denominator)
    return value

if __name__ == "__main__":
    # The problem asks for the exact value of l(d)
    # The derived formula is l(d) = ln((sqrt(d) - 1) / (sqrt(d) + 1))
    
    # Let's present the formula as the answer.
    # The equation involves variables and functions, so we will print it as a string.
    # Using 'd' as a symbolic representation.
    
    print("The exact value of l(d) is given by the equation:")
    print("l(d) = ln((sqrt(d) - 1) / (sqrt(d) + 1))")
    
    # As an example, let's calculate and print the value for a specific d, e.g., d=4
    d_example = 4
    l_value = calculate_l(d_example)
    sqrt_d_val = math.sqrt(d_example)
    
    print(f"\nFor example, when d = {d_example}:")
    print(f"l({d_example}) = ln((sqrt({d_example}) - 1) / (sqrt({d_example}) + 1))")
    print(f"l({d_example}) = ln(({sqrt_d_val:.2f} - 1) / ({sqrt_d_val:.2f} + 1))")
    print(f"l({d_example}) = ln(({sqrt_d_val-1}) / ({sqrt_d_val+1}))")
    print(f"l({d_example}) = ln(1/3)")
    print(f"l({d_example}) = {l_value}")
