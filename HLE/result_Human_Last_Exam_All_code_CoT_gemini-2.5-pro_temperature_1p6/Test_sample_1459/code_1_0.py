import math

def calculate_gh_distance_interval_circle():
    """
    Calculates the Gromov-Hausdorff distance between the interval [0,1]
    and the unit circle with the intrinsic metric.
    """
    # Define the parameters of the two metric spaces
    # 1. The interval [0, 1]
    L = 1.0  # Length of the interval

    # 2. The unit circle
    R = 1.0  # Radius of the circle

    # The Gromov-Hausdorff distance between an interval of length L
    # and a circle of radius R is given by the formula:
    # d_GH = (1/2) * sqrt(L^2 + (pi*R)^2)

    # Calculate the components of the formula step-by-step
    L_squared = L**2
    pi_R = math.pi * R
    pi_R_squared = pi_R**2
    sum_of_squares = L_squared + pi_R_squared
    sqrt_of_sum = math.sqrt(sum_of_squares)
    gh_distance = 0.5 * sqrt_of_sum

    # Print the explanation and the step-by-step calculation
    print("The Gromov-Hausdorff distance between the interval [0, 1] (length L=1) and the unit circle (radius R=1) is calculated using the formula:")
    print("d_GH = (1/2) * sqrt(L^2 + (pi*R)^2)\n")
    print("Substituting L = 1.0 and R = 1.0:")
    
    # Printing each number in the final equation
    print(f"d_GH = (1/2) * sqrt({L:.1f}^2 + (pi * {R:.1f})^2)")
    print(f"d_GH = (1/2) * sqrt({L_squared:.4f} + ({math.pi:.4f} * {R:.1f})^2)")
    print(f"d_GH = (1/2) * sqrt({L_squared:.4f} + {pi_R_squared:.4f})")
    print(f"d_GH = (1/2) * sqrt({sum_of_squares:.4f})")
    print(f"d_GH = 0.5 * {sqrt_of_sum:.4f}")
    print(f"d_GH = {gh_distance:.4f}\n")
    print(f"The final calculated distance is: {gh_distance}")

calculate_gh_distance_interval_circle()