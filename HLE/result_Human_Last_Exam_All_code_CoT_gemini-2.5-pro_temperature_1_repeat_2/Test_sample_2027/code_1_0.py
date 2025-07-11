import math

def calculate_l_d(d):
    """
    This function calculates the exact value of l(d) based on the derived formula.

    Args:
        d (int): The dimension of the Poincaré disk, must be an integer >= 2.

    Returns:
        float: The exact value of l(d).
    """
    if not isinstance(d, int) or d < 2:
        raise ValueError("d must be an integer greater than or equal to 2.")

    # The exact value is given by the formula: l(d) = ln((sqrt(d) - 1) / (sqrt(d) + 1))
    sqrt_d = math.sqrt(d)
    
    # These are the numbers in the final simplified equation's fraction
    number_in_numerator = 1
    number_in_denominator = 1

    numerator_val = sqrt_d - number_in_numerator
    denominator_val = sqrt_d + number_in_denominator
    
    result = math.log(numerator_val / denominator_val)
    
    print(f"The final equation for l(d) is: l(d) = ln((sqrt(d) - {number_in_numerator}) / (sqrt(d) + {number_in_denominator}))")
    print(f"For the specific case where d = {d}:")
    print(f"l({d}) = ln((sqrt({d}) - {number_in_numerator}) / (sqrt({d}) + {number_in_denominator}))")
    print(f"l({d}) = ln({numerator_val} / {denominator_val})")
    print(f"The exact value is: {result}")


# --- Main execution ---
# You can change the value of d here.
d_value = 4
try:
    calculate_l_d(d_value)
    # The final answer is the derived expression.
    # For d=4, the value is ln((2-1)/(2+1)) = ln(1/3) approx -1.0986
    final_answer = math.log((math.sqrt(d_value) - 1) / (math.sqrt(d_value) + 1))
    # We are asked for the exact value of l(d), which is the symbolic expression.
    # Here is an example numeric result for d=4
    print(f"<<<l({d_value}) = {final_answer}>>>")
except ValueError as e:
    print(e)
