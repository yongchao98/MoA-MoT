import math

def solve_for_cp():
    """
    Calculates the length of CP given the length of AB.
    """
    # The length of segment AB is given.
    ab_length = 18468

    # The derived formula is: CP = AB * (sqrt(5) + 1) / 2
    # Let's define the numbers in the equation.
    num_ab = ab_length
    num_sqrt = 5
    num_one = 1
    num_two = 2

    # Calculate the final result.
    cp_length = num_ab * (math.sqrt(num_sqrt) + num_one) / num_two

    # Print the equation with the numbers plugged in, as requested.
    print("Based on the geometric derivation, the formula for the length of CP is:")
    print("CP = AB * (sqrt(5) + 1) / 2")
    print("\nPlugging in the given value for AB:")
    print(f"CP = {num_ab} * (sqrt({num_sqrt}) + {num_one}) / {num_two}")

    # Print the final calculated length, rounded to 4 decimal points.
    print(f"\nThe calculated length of CP is:")
    print(f"{cp_length:.4f}")

solve_for_cp()