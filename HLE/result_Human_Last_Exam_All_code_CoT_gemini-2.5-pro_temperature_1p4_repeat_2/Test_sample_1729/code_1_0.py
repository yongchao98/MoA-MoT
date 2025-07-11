import math

def print_probability_formula():
    """
    This function prints the derived formula for the probability P_m.
    """
    
    # The number of valid pairs (i, j) is C(m+2, 2).
    # Its formula is (m+1)(m+2)/2.
    numerator_formula = "(m + 1)(m + 2)"
    
    # The total number of pairs is C(4m+2, 2).
    # Its formula is (2m+1)(4m+1).
    # The probability is the ratio of valid pairs to total pairs.
    # P_m = [ (m+1)(m+2)/2 ] / [ (2m+1)(4m+1) ]
    denominator_formula = "2(2m + 1)(4m + 1)"

    print("The probability P_m is given by the formula:")
    print(f"      {numerator_formula}")
    print(f"P_m = —————————————————————")
    print(f"      {denominator_formula}")
    
    # The prompt requests to output each number in the final equation.
    # The equation is symbolic, but we can explain its components.
    print("\nWhere the formula components are:")
    print(f"Numerator, number of valid pairs = C(m+2, 2) = {numerator_formula} / 2")
    print(f"Denominator, total number of pairs = C(4m+2, 2) = (2m+1)(4m+1)")

# Execute the function to display the result.
print_probability_formula()