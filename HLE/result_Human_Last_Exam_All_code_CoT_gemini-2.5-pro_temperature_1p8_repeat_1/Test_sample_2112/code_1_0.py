import math

def solve_for_r0():
    """
    Calculates the radial distance r_0 > 15 where the gravitational potential vanishes.
    The derivation shows that r_0 must satisfy the equation:
    17 * r_0 = 226 + 49 * sqrt(2)
    """

    # The equation for r_0 is of the form: c1 * r_0 = c2 + c3 * sqrt(2)
    c1 = 17
    c2 = 226
    c3 = 49

    print("The problem reduces to solving a specific algebraic equation for r_0.")
    print("The final derived equation for r_0 is:")
    # We output the numbers in the final equation.
    print(f"{c1} * r_0 = {c2} + {c3} * sqrt(2)")
    print("")

    # Calculate the value of r_0
    sqrt_2 = math.sqrt(2)
    r_0 = (c2 + c3 * sqrt_2) / c1

    print(f"The calculated value for the radial distance r_0 is: {r_0}")

if __name__ == "__main__":
    solve_for_r0()