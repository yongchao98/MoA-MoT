import math

def solve_dna_probability():
    """
    Calculates the limiting probability P(n) as n -> infinity.

    The problem states that replication fails if for any base x (from 0 to 7),
    the count of that base, N_x, satisfies N_x mod 8 = x.
    The probability of success is the probability that N_x mod 8 != x for all x.

    As n -> infinity, the vector of counts (N_0 mod 8, ..., N_7 mod 8) becomes
    uniformly distributed over the 8^8 possible states.

    We need to count the number of "success" states, where for each position x,
    the value is not equal to x. For each of the 8 positions, there are 7 valid values.
    So, the number of success states is 7^8.

    The total number of states is 8^8.
    The limiting probability is therefore (7^8) / (8^8) = (7/8)^8.
    """
    
    numerator_base = 7
    denominator_base = 8
    exponent = 8

    # Calculate the numerator and the denominator of the probability
    numerator = numerator_base ** exponent
    denominator = denominator_base ** exponent

    # Calculate the final probability as a float
    probability = numerator / denominator

    # Print the final equation with all its components as requested
    print(f"The limiting probability P(n) as n -> infinity is given by the expression: ({numerator_base}/{denominator_base})^{exponent}")
    print(f"P(n->inf) = {numerator_base}^{exponent} / {denominator_base}^{exponent} = {numerator} / {denominator} = {probability}")

solve_dna_probability()