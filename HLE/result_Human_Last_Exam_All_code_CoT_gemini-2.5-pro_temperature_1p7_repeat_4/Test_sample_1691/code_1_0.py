import math

def solve():
    """
    This function calculates the analytical approximation for the integral I(epsilon)
    for the small epsilon regime.
    """
    
    # Parameters from the dominant term in the denominator, 9.0 * x^5.0
    c = 9.0
    p = 5.0
    
    # The analytical approximation has the form I(epsilon) = C * epsilon^a.
    # We calculate the exponent 'a' and the coefficient 'C' based on the general formula.
    # a = (1/p) - 1
    # C = (pi / (p * sin(pi/p))) * c^(-1/p)
    
    # Calculate the exponent 'a'
    a = (1.0 / p) - 1.0
    
    # Calculate the coefficient 'C'
    # The calculation is split into parts for clarity.
    c_inv_p = c**(-1.0/p)
    sine_term = math.sin(math.pi / p)
    
    C = (math.pi / (p * sine_term)) * c_inv_p
    
    print("The analytical formula for I(epsilon) for small epsilon is of the form:")
    print("I(epsilon) = C * epsilon^a")
    print("\nBased on the integral, the parameters of the dominant term (c * x^p) are:")
    print(f"c = {c}")
    print(f"p = {p}")
    print("\nThe calculated values for the constants in the formula are:")
    print(f"Exponent a = {a}")
    print(f"Coefficient C = {C}")
    print("\nThus, the final approximate analytical formula is:")
    print(f"I(epsilon) = {C} * epsilon^({a})")

solve()