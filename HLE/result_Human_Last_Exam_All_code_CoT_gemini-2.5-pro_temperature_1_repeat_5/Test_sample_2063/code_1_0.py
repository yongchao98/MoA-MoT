import math

def solve_quantum_dna_probability():
    """
    Calculates the limiting probability P(n) as n -> infinity.

    The problem states that replication succeeds if for all 8 bases (x=0 to 7),
    the count of that base (c_x) does not satisfy c_x mod 8 = x.

    As n approaches infinity, the probability that any single one of these
    conditions is met (c_x mod 8 = x) approaches 1/8. The probability that it is not
    met is therefore 7/8.

    While the counts (c_x) are not truly independent because their sum must be n,
    the limiting probability, when averaged over all possibilities of n mod 8,
    converges to the value as if they were independent.

    Therefore, the probability that *all 8* success conditions are met is the
    product of their individual probabilities.
    """

    base = 7
    denominator = 8
    exponent = 8

    # Calculate the result
    result = (base / denominator) ** exponent

    # Print the equation and the final value
    print(f"The limiting probability is given by the expression ({base}/{denominator})^{exponent}.")
    print(f"Calculating the values: {base}^{exponent} = {base**exponent} and {denominator}^{exponent} = {denominator**exponent}.")
    print(f"The final equation is: {base**exponent} / {denominator**exponent} = {result}")

solve_quantum_dna_probability()