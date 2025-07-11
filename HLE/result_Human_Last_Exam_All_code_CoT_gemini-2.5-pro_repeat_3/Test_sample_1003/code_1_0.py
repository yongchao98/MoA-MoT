import math

def solve_star_angle_ratio():
    """
    This function calculates the value of (1 - cos(theta_14')) / (1 - cos(theta_34')).

    The solution is based on the Lorentz invariance of the 4-vector dot product, which leads
    to the relation D_i * D_j * (1 - cos(theta'_ij)) = constant for any pair of stars i, j.
    From the problem statement, we deduce that the desired ratio simplifies to D_3 / D_1, where
    D_i is the Doppler factor for star i.

    Further analysis of the given angles in the second frame leads to the equation:
    D_1 = D_3 * (1 + 1/sqrt(2)).
    Therefore, the ratio D_3 / D_1 is 1 / (1 + 1/sqrt(2)).

    This simplifies to 2 - sqrt(2).
    """

    # The numbers in the final equation
    num1 = 2
    sqrt_of_2 = math.sqrt(2)

    # Calculate the final value
    final_value = num1 - sqrt_of_2

    print("The final expression for the ratio is 2 - sqrt(2).")
    print("Let's calculate the numerical value:")
    # Using f-string to format the output and show the numbers in the equation
    print(f"{num1} - {sqrt_of_2:.5f} = {final_value:.5f}")
    print(f"\nThe exact value is 2 - sqrt(2).")
    print(f"The numerical value is approximately {final_value}")

solve_star_angle_ratio()