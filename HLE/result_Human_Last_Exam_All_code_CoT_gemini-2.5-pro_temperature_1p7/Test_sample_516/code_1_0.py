import math

def solve_expected_value():
    """
    This function calculates the expected value E[T] based on the derived formula.
    The formula is E[T] = 3 - 2 * e^(1/3).
    """

    # The numbers in the final equation E[T] = a - b * e^(c/d)
    a = 3
    b = 2
    c = 1
    d = 3

    # Calculate the result
    result = a - b * math.exp(c / d)

    # Output the components of the final equation and the result
    print(f"The derived equation for the expected value E[T] is of the form: a - b * e^(c/d)")
    print(f"The numbers in the equation are:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    print(f"d = {d}")
    print(f"\nThe final equation is: {a} - {b} * e^({c}/{d})")
    print(f"\nThe calculated expected value E[T] is: {result}")

if __name__ == "__main__":
    solve_expected_value()
