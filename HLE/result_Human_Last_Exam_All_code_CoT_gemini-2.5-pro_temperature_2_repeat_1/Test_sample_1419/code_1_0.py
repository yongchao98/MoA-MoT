def solve_fixed_point_coupling():
    """
    Calculates and prints the leading order expression for the Wilson-Fisher
    fixed point coupling u* in phi^4 theory near d=4 dimensions.
    """

    # The problem setup
    print("In phi^4 theory in d = 4 - ε dimensions, the one-loop beta function for the coupling u is:")
    # The numbers in this equation are 1 (implicit in εu), and 3.
    print("β(u) = -1 * ε * u + 3 * u^2")
    print("-" * 40)

    # The fixed-point condition
    print("Fixed points u* are found by solving β(u*) = 0:")
    print("-ε*u* + 3*(u*)^2 = 0")
    print("\nFactoring out u* gives:")
    print("u* * (-ε + 3*u*) = 0")
    print("-" * 40)


    # The solutions
    print("This equation has two solutions:")
    print("1. The Gaussian fixed point: u* = 0")
    print("2. The non-trivial Wilson-Fisher fixed point, which is found by solving:")
    print("-ε + 3*u* = 0")
    print("\nRearranging the terms to solve for u*:")
    print("3*u* = ε")
    print("-" * 40)

    # The final answer
    # The final equation is u* = ε / 3. The numbers involved are 1 (implicit for u* and ε) and 3.
    # The following print statement outputs the full equation including all numbers.
    print("The leading order expression for the fixed point coupling u* is:")
    print("u* = (1 * ε) / 3")


if __name__ == '__main__':
    solve_fixed_point_coupling()
