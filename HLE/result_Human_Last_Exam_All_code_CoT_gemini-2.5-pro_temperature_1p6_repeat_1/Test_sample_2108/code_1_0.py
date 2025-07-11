import math

def calculate_max_ratio():
    """
    Calculates the maximum achievable ratio of bidirectional conical power to line intensity.
    
    The final expression for the ratio is derived from the physics of electromagnetic
    radiation from a specific rotating quadrupole charge distribution. The optimal
    configuration involves rotating the sphere about a non-null principal axis
    of its quadrupole tensor.
    
    The ratio is given by the analytical expression:
    Ratio = 4 * pi * (28/15 - 83 / (60 * sqrt(2)))
    """
    
    # Numbers used in the final equation
    num_4 = 4
    pi_val = math.pi
    num_28 = 28
    num_15 = 15
    num_83 = 83
    num_60 = 60
    num_2 = 2
    sqrt_2 = math.sqrt(num_2)

    print("The final calculation is based on the equation: R = A * pi * (B/C - D / (E * sqrt(F)))")
    print(f"Where the numbers in the equation are:")
    print(f"A = {num_4}")
    print(f"pi = {pi_val}")
    print(f"B = {num_28}")
    print(f"C = {num_15}")
    print(f"D = {num_83}")
    print(f"E = {num_60}")
    print(f"F = {num_2}")
    print("-" * 20)

    # Perform the calculation
    term1 = num_28 / num_15
    term2 = num_83 / (num_60 * sqrt_2)
    
    result = num_4 * pi_val * (term1 - term2)
    
    print(f"The calculated maximum ratio is: {result}")

calculate_max_ratio()