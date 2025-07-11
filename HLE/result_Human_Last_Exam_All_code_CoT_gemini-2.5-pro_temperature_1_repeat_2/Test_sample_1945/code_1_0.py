import math

def calculate_limit_probability():
    """
    Calculates and prints the exact value of lim Pr(n) as n -> infinity.
    The limit is derived analytically to be 3 / pi^2.
    This code demonstrates the calculation of that final value.
    """
    # The final equation for the limit is L = numerator / denominator.
    
    # From the analytical derivation, the numerator is 3.
    numerator = 3
    
    # The denominator from the derivation is pi squared.
    pi_value = math.pi
    denominator = pi_value**2
    
    # Calculate the final result
    limit_value = numerator / denominator
    
    print("The analytical derivation shows the limit of the probability is given by the equation: 3 / (pi^2)")
    print("\n--- Calculation Breakdown ---")
    print(f"Numerator: {numerator}")
    print(f"Value of pi: {pi_value}")
    print(f"Denominator (pi^2): {denominator}")
    print("--------------------------")
    print(f"The exact value of the limit is {limit_value}")

calculate_limit_probability()
