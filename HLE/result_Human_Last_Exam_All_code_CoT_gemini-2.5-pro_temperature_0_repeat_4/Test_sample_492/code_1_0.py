import math

def solve_for_c():
    """
    This function calculates the exact value of the constant c, which represents
    the time of emergence of the giant connected component in the described model.

    The derivation shows that the average degree k(t) of the graph at time t is given by:
    k(t) = t^2 / 3

    The giant component emerges when the average degree k(c) is equal to 1.
    This gives the critical equation: c^2 / 3 = 1.
    """

    # The critical equation is c^2 / denominator = numerator
    numerator = 1
    denominator = 3

    # We solve for c^2
    # c^2 = numerator * denominator
    c_squared = numerator * denominator

    # The value of c is the square root of c_squared
    c = math.sqrt(c_squared)

    print("The critical equation for the emergence of the giant component is:")
    print(f"c^2 / {denominator} = {numerator}")
    print("\nSolving for c:")
    print(f"c^2 = {numerator} * {denominator}")
    print(f"c^2 = {c_squared}")
    print(f"c = sqrt({c_squared})")
    print("\nThe exact value of c is the square root of 3.")
    print(f"c = {c}")

if __name__ == "__main__":
    solve_for_c()