def solve_hypercube_problem():
    """
    Finds the values of n for which it is possible to have only one gift left.
    The condition for solvability is that n must be 1 modulo 7.
    """
    
    # The range for n is 7 <= n <= 55.
    min_n = 7
    max_n = 55
    
    # The length of the move is 7.
    move_length = 7
    
    # According to the mathematical derivation, a solution exists if and only if
    # n % move_length == 1.
    possible_n_values = []
    for n in range(min_n, max_n + 1):
        if n % move_length == 1:
            possible_n_values.append(n)
            
    # The problem asks for the output to be printed.
    # The "final equation" part of the prompt is interpreted as printing the final result set.
    # We will print the numbers separated by spaces.
    print(*possible_n_values)

solve_hypercube_problem()