import math

def solve_lewis_weight_constant():
    """
    Calculates the smallest constant c for the inequality involving Lewis weights.
    """
    
    # We fix the dimension d and parameter p for demonstration purposes.
    # The user can change these values to see the result for other settings.
    d = 10
    p = 1.5

    print(f"For a fixed dimension d = {d} and parameter p = {p}:")
    
    # The smallest constant c is given by the formula c = d**|1/2 - 1/p|.
    
    # Calculate the exponent
    exponent_val = 0.5 - (1/p)
    abs_exponent_val = abs(exponent_val)

    # Calculate the constant c
    c = d ** abs_exponent_val
    
    print("\nThe smallest constant c is determined by the formula: c = d**|1/2 - 1/p|")
    print("\nCalculation steps:")
    # Using f-string formatting to show each number in the equation
    print(f"c = {d} ** |1/2 - 1/{p}|")
    print(f"c = {d} ** |{0.5} - {1/p}|")
    print(f"c = {d} ** |{exponent_val}|")
    print(f"c = {d} ** {abs_exponent_val}")
    print(f"Final value of c: {c}")

solve_lewis_weight_constant()