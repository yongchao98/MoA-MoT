import math

def main():
    """
    This script derives and calculates the area ratio for the described geometric construction.
    """
    # The general formula for the ratio of the areas, derived from the plan.
    # Ratio = cos^2(pi / (2*n)) / cos(pi / n)
    
    print("The general formula for how many times larger the area of the n-sided polygon is compared to the 2n-sided polygon is:")
    print("Ratio = cos^2(pi / (2*n)) / cos(pi / n)")
    
    print("\n--------------------------------------------------")
    print("Let's verify this with the example from the problem: a hexagon (2n=6) inside a triangle (n=3).")
    
    n = 3
    
    # Angles in radians for the formula
    angle_in_cos_squared = math.pi / (2 * n)
    angle_in_cos = math.pi / n
    
    # Calculate the values of the trigonometric functions
    cos_val_for_squared = math.cos(angle_in_cos_squared)
    cos_val = math.cos(angle_in_cos)
    
    # Calculate the final ratio
    ratio = (cos_val_for_squared ** 2) / cos_val
    
    print(f"\nFor n = {n}:")
    print(f"The equation is: Ratio = cos^2(pi / (2 * {n})) / cos(pi / {n})")
    print(f"1. Calculate the term in the numerator: cos(pi / {2*n}) = cos({angle_in_cos_squared:.4f} rad) = {cos_val_for_squared:.4f}")
    print(f"2. Square the result: ({cos_val_for_squared:.4f})^2 = {cos_val_for_squared**2:.4f}")
    print(f"3. Calculate the term in the denominator: cos(pi / {n}) = cos({angle_in_cos:.4f} rad) = {cos_val:.4f}")
    print(f"4. Divide the numerator by the denominator: {cos_val_for_squared**2:.4f} / {cos_val:.4f} = {ratio:.4f}")
    
    print(f"\nThe final ratio is {ratio:.1f}, which is exactly 3/2, matching the example.")

if __name__ == "__main__":
    main()