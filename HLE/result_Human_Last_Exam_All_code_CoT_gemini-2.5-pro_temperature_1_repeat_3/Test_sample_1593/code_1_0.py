def solve_complexity():
    """
    This function calculates the complexity tuple (a, b, c) based on the analytical derivation.

    The query complexity Q for both regimes was determined to be Theta(N * log(N)).
    We are asked to represent this in the form sqrt(N^a * (log N)^b * (log log N)^c).

    Equating the two forms:
    N * log(N) = (N^a * (log N)^b * (log log N)^c)^0.5

    Squaring both sides:
    (N * log(N))^2 = N^a * (log N)^b * (log log N)^c
    N^2 * (log N)^2 * (log log N)^0 = N^a * (log N)^b * (log log N)^c

    By comparing the exponents for N, log N, and log log N, we find the values for a, b, and c.
    """

    # For the N term: N^2 = N^a  => a = 2
    a = 2
    # For the log(N) term: (log N)^2 = (log N)^b => b = 2
    b = 2
    # For the log(log(N)) term: (log log N)^0 = (log log N)^c => c = 0
    c = 0

    # The result (a,b,c) is the same for both regimes.
    regime1_tuple = (a, b, c)
    regime2_tuple = (a, b, c)

    # Print the final answer in the specified format "(a,b,c),(a,b,c)".
    # The requirement to "output each number in the final equation" is interpreted
    # as programmatically constructing the final output string from the derived numbers.
    print(f"({regime1_tuple[0]},{regime1_tuple[1]},{regime1_tuple[2]}),({regime2_tuple[0]},{regime2_tuple[1]},{regime2_tuple[2]})")

solve_complexity()