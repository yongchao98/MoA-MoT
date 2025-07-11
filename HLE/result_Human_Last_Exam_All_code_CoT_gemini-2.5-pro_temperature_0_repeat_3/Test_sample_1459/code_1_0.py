import math

def calculate_gromov_hausdorff_distance():
    """
    Calculates the Gromov-Hausdorff distance between the interval [0,1]
    and the unit circle with the intrinsic metric.
    """
    # Length of the interval [0,1]
    l = 1.0

    # Circumference of the unit circle (radius=1)
    L = 2 * math.pi

    # The formula for the distance is d_GH = (1/4) * sqrt((L - 2*l)^2 + l^2)
    # We first calculate the term inside the square root
    term_inside_sqrt = (L - 2 * l)**2 + l**2

    # Then we calculate the final distance
    distance = (1 / 4) * math.sqrt(term_inside_sqrt)

    # As requested, we print the final equation with each number.
    # The numbers in the equation are 1, 4, L, 2, and l.
    print("The Gromov-Hausdorff distance is calculated using the formula: d_GH = (1/4) * sqrt((L - 2*l)^2 + l^2)")
    print("\nSubstituting L (circumference of unit circle) and l (length of interval):")
    print(f"d_GH = (1 / 4) * sqrt(( {L} - 2 * {l} )^2 + {l}^2)")

    # Print the final numerical result
    print("\nThe calculated distance is:")
    print(distance)

if __name__ == "__main__":
    calculate_gromov_hausdorff_distance()