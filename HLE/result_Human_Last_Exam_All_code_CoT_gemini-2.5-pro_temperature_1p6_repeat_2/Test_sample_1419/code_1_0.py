def solve_fixed_point_coupling():
    """
    Derives and prints the leading order expression for the Wilson-Fisher
    fixed point coupling u* in ϕ⁴ theory near four dimensions.
    """
    # Step 1: State the one-loop beta function
    print("Step 1: The starting point is the one-loop beta function for the dimensionless coupling 'u'")
    print("in ϕ⁴ theory, dimensionally regularized in d = 4 - ε dimensions.")
    print("The beta function describes how the coupling changes with energy scale.")
    print("To leading order in u and ε, it is given by:")
    print("β(u) = -ε·u + (3 / (16·π²))·u²")
    print("-" * 50)

    # Step 2: Define the fixed point condition
    print("Step 2: A fixed point, denoted u*, is a point where the coupling ceases to change,")
    print("meaning the beta function is zero. This corresponds to a scale-invariant theory.")
    print("So, we set β(u*) = 0:")
    print("-ε·u* + (3 / (16·π²))·(u*)² = 0")
    print("-" * 50)

    # Step 3: Solve the equation for the non-trivial fixed point
    print("Step 3: We solve this equation for u*.")
    print("The equation can be rewritten by factoring out u*:")
    print("u* · [-ε + (3 / (16·π²))·u*] = 0")
    print("\nThis gives two possible solutions:")
    print("1. The trivial (Gaussian) fixed point: u* = 0")
    print("2. The non-trivial (Wilson-Fisher) fixed point, found by setting the term in brackets to zero.")
    print("\nSolving for the non-trivial fixed point:")
    print("-ε + (3 / (16·π²))·u* = 0")
    print("(3 / (16·π²))·u* = ε")
    print("3·u* = 16·π²·ε")
    print("u* = (16·π² / 3)·ε")
    print("-" * 50)

    # Step 4: Final answer
    print("Step 4: The final leading order expression for the Wilson-Fisher fixed point coupling is:")
    numerator = 16
    denominator = 3
    # As requested, the final equation with each number explicitly shown:
    print(f"u* = ({numerator}·π² / {denominator})·ε")

if __name__ == "__main__":
    solve_fixed_point_coupling()