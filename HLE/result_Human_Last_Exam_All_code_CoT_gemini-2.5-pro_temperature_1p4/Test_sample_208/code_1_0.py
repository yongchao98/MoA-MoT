import math

def calculate_packing_radius():
    """
    Calculates the radius of a circle that tightly packs 14 unit circles.
    """
    # The radius of the small circles
    r = 1.0

    # The best-known packing for 14 circles gives the radius of the
    # containing circle R with the formula: R = r * (1 + sqrt(2) + sqrt(6))
    
    # Calculate the necessary square roots
    val_sqrt_2 = math.sqrt(2)
    val_sqrt_6 = math.sqrt(6)
    
    # Calculate the final radius R
    R = r * (1 + val_sqrt_2 + val_sqrt_6)

    # As requested, printing the final equation with each number.
    print("The formula for the radius 'R' of the large circle is:")
    print("R = r * (1 + sqrt(2) + sqrt(6))")
    print("\nSubstituting the values:")
    print(f"r = {r}")
    print(f"sqrt(2) = {val_sqrt_2}")
    print(f"sqrt(6) = {val_sqrt_6}")
    print("\nThe equation becomes:")
    print(f"R = {r} * (1 + {val_sqrt_2} + {val_sqrt_6})")
    
    # Calculate and print the final result before rounding
    total_ratio = 1 + val_sqrt_2 + val_sqrt_6
    print(f"R = {r} * ({total_ratio})")
    print(f"R = {R}")

    # Format the result to 4 significant digits.
    # The value is 4.8637... The 5th significant digit is 7, so we round up.
    # The result becomes 4.864.
    final_answer_str = f"{R:.4g}"

    print(f"\nThe radius of the large circle up to 4 significant digits is {final_answer_str}.")

calculate_packing_radius()