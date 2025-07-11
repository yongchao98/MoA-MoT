def solve_limiting_probability():
    """
    Calculates the limiting probability P(n) as n -> infinity.

    The probability P(n) is for the successful replication of a DNA sequence
    of length n, where success is defined as the count of each of the 8 bases `x`
    (c_x) not being congruent to `x` modulo 8.

    As n approaches infinity, the events {c_x mod 8 = x} become
    independent, and the probability of each such event approaches 1/8.
    The probability of the complementary event {c_x mod 8 != x} is 1 - 1/8 = 7/8.

    The overall probability of success is the product of these probabilities for
    all 8 bases, which is (7/8)^8.
    """

    # Define the parameters for the final equation
    numerator = 7
    denominator = 8
    exponent = 8

    # Calculate the components of the final answer
    result_numerator = numerator ** exponent
    result_denominator = denominator ** exponent
    decimal_result = result_numerator / result_denominator

    # Output the explanation and the final equation step-by-step
    print("The closed-form expression for the limiting probability P is (7/8)^8.")
    print("This is derived from the product of probabilities for each of the 8 bases not meeting the collapse condition.")
    print("For each base x, P(c_x mod 8 != x) -> 1 - 1/8 = 7/8.")
    print("\nThe calculation is as follows:")

    # Output each number in the final equation as requested
    print(f"P = ({numerator}/{denominator})^{exponent}")
    print(f"P = {numerator}^{exponent} / {denominator}^{exponent}")
    print(f"P = {result_numerator} / {result_denominator}")
    print(f"P ≈ {decimal_result}")

solve_limiting_probability()