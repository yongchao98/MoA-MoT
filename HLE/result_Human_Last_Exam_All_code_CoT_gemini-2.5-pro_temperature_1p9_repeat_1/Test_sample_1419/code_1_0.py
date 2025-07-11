import math

def get_fixed_point_expression():
    """
    This function derives and prints the leading order expression for the
    fixed point coupling u* in phi^4 theory near d=4 dimensions.
    """
    
    # In d = 4 - ε dimensions, the one-loop beta function is:
    # β(u) = -εu + (3u²) / (16π²)
    
    # A fixed point u* is found by solving β(u*) = 0.
    # We are interested in the non-trivial (Wilson-Fisher) fixed point.
    # -εu* + (3(u*)²) / (16π²) = 0
    #  => u* = (16π²ε) / 3
    
    # Constants in the expression
    numerator = 16
    denominator = 3
    pi_symbol = "π"
    epsilon_symbol = "ε"
    
    print("The leading order expression for the fixed point coupling u* is derived from the beta function β(u).")
    print("The one-loop beta function in d = 4 - ε dimensions is: β(u) = -εu + (3u²)/(16π²)")
    print("\nThe non-trivial fixed point u* (Wilson-Fisher fixed point) is found by solving β(u*) = 0:")
    
    # Output the final equation, showing each number.
    print(f"u* = ({numerator} * {pi_symbol}² / {denominator}) * {epsilon_symbol}")

if __name__ == "__main__":
    get_fixed_point_expression()