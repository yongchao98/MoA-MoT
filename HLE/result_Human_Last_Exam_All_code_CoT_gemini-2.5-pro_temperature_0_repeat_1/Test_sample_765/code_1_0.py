import math
import random

def is_in_R(a, b):
    """
    Checks if a point z = a + bi is in the region R.
    The condition is floor(|z|) = |floor(z)|.
    """
    # Let na = floor(a) and nb = floor(b)
    na = math.floor(a)
    nb = math.floor(b)

    # Calculate the squared magnitude of the floored complex number
    k_squared = na*na + nb*nb
    
    # The magnitude |floor(z)| must be an integer for the main equation to hold,
    # since the left side floor(|z|) is always an integer.
    # So, k = sqrt(k_squared) must be an integer.
    k = math.sqrt(k_squared)
    if k != math.floor(k):
        return False

    # Now check the full condition
    k_int = int(k)
    magnitude_z = math.sqrt(a*a + b*b)
    floor_magnitude_z = math.floor(magnitude_z)

    return floor_magnitude_z == k_int

def estimate_area():
    """
    Estimates the area of R using Monte Carlo simulation.
    """
    num_points = 10_000_000  # Total number of random points to generate
    count_in_R = 0          # Counter for points inside R

    # The square is defined by 0 <= a <= 6 and 0 <= b <= 6
    for _ in range(num_points):
        a = random.uniform(0, 6)
        b = random.uniform(0, 6)
        if is_in_R(a, b):
            count_in_R += 1
    
    total_domain_area = 6 * 6
    estimated_area = (count_in_R / num_points) * total_domain_area
    
    print("The area of R is estimated using the Monte Carlo method.")
    print("Area = (Points in R / Total Points) * Total Area")
    print(f"Area = ({count_in_R} / {num_points}) * {total_domain_area}")
    print(f"Calculated Area: {estimated_area:.2f}")

if __name__ == '__main__':
    estimate_area()
