def solve_hypercube_puzzle():
    """
    Finds and prints the integer values for n in the range [7, 55] for which it is
    possible to leave a single gift on a 5D hypercube of side n.

    The mathematical analysis of this puzzle using linear algebra over the finite
    field F_2 shows that a final state with one gift is reachable from an initial
    state with all gifts if and only if n is congruent to 1 modulo 7.
    """

    # The problem specifies the range for n.
    start_n = 7
    end_n = 55
    
    # According to the derivation, n must be of the form 7k + 1.
    divisor = 7
    required_remainder = 1
    
    valid_n_values = []
    for n in range(start_n, end_n + 1):
        if n % divisor == required_remainder:
            valid_n_values.append(n)
            
    # The problem asks for the values in increasing order.
    # The loop generates them in increasing order, so no sorting is needed.
    print("The values of n for which it is possible to reach the state with only one gift are:")
    
    # We print each number, separated by commas, as requested.
    output_str = ", ".join(map(str, valid_n_values))
    print(output_str)

solve_hypercube_puzzle()