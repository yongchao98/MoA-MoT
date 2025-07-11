import math

def solve_complexity():
    """
    This function determines the query complexity for the two specified regimes
    and prints the result in the (a,b,c) format.
    The logic is based on the theoretical derivation explained above.
    """

    # For both regimes, the analysis shows the optimal sorting strategy has a
    # query complexity of Theta(N * log(N)).

    # We need to express this complexity in the form:
    # Theta(sqrt(N^a * (logN)^b * (loglogN)^c))

    # Setting the two forms equal:
    # sqrt(N^a * (logN)^b * (loglogN)^c) = N * logN

    # Squaring both sides of the equation:
    # N^a * (logN)^b * (loglogN)^c = (N * logN)^2
    # N^a * (logN)^b * (loglogN)^c = N^2 * (logN)^2

    # By comparing the exponents of N, logN, and loglogN, we can find the
    # values for a, b, and c.
    # For N: a = 2
    # For logN: b = 2
    # For loglogN: c = 0
    a = 2
    b = 2
    c = 0

    # The (a,b,c) tuple is the same for both regimes.
    regime1_abc = (a, b, c)
    regime2_abc = (a, b, c)

    # We print the numbers that form the final answer tuple, as requested.
    # The final output is a string containing two tuples for the two regimes.
    print(f"({regime1_abc[0]},{regime1_abc[1]},{regime1_abc[2]}),({regime2_abc[0]},{regime2_abc[1]},{regime2_abc[2]})")

solve_complexity()