import math

def solve_radiation_ratio():
    """
    This script calculates the maximum achievable ratio of bidirectional conical
    power to line intensity for a rotating charged sphere as described.

    The derivation shows that the problem reduces to analyzing the radiation
    pattern I(theta) = C * (1 + cos^2(theta)), where theta is the angle from
    the rotation axis. The maximum ratio is found by aligning the cone with
    the rotation axis (to maximize power) and choosing the line perpendicular
    to the rotation axis (to minimize intensity).

    The final ratio is given by the formula: (pi/3) * (16 - 7*sqrt(2)).
    This script calculates and prints the value of this expression.
    """
    print("This program calculates the maximum ratio of bidirectional conical power to line intensity.")
    print("The derived formula for the ratio is: (pi/3) * (16 - 7 * sqrt(2))")
    print("\n--- Breakdown of the calculation ---")

    # The final equation is Ratio = (pi/3) * (16 - 7 * sqrt(2))
    # We will output each number/component of this equation.

    pi = math.pi
    print(f"Value of pi: {pi}")

    sqrt_2 = math.sqrt(2)
    print(f"Value of sqrt(2): {sqrt_2}")

    term1 = 7 * sqrt_2
    print(f"Value of '7 * sqrt(2)': {term1}")

    term2 = 16 - term1
    print(f"Value of '(16 - 7 * sqrt(2))': {term2}")

    term3 = pi / 3
    print(f"Value of 'pi / 3': {term3}")

    final_ratio = term3 * term2
    print("\n--- Final Result ---")
    print(f"The maximum achievable ratio is (pi/3) * (16 - 7*sqrt(2))")
    print(f"Maximum Ratio = {final_ratio}")

if __name__ == "__main__":
    solve_radiation_ratio()
    # The final numerical result for the ratio.
    # The calculation is (math.pi/3) * (16 - 7 * math.sqrt(2))
    final_answer = (math.pi/3) * (16 - 7 * math.sqrt(2))
    # To conform to the output format, we put the answer here.
    # print(f"<<<{final_answer}>>>")