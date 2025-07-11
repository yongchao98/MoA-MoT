def find_solvable_hypercube_sizes():
    """
    Calculates and prints the values of n for which the hypercube puzzle is solvable.

    The state of the n x n x n x n x n hypercube is represented by a vector
    over the field F_2. A move consists of flipping the states of 7
    consecutive positions along an axis. The puzzle starts with all gifts present (all 1s)
    and the goal is to reach a state with a single gift (a single 1).

    The problem is solvable if and only if n is not divisible by 7.
    This code iterates through the given range for n and prints the values
    that satisfy this condition.
    """
    
    min_n = 7
    max_n = 55
    
    solvable_n_values = []
    
    for n in range(min_n, max_n + 1):
        # The condition for solvability is that n is not divisible by 7.
        if n % 7 != 0:
            solvable_n_values.append(n)
            
    print("The values of n for which it is possible to have only one gift are:")
    # The final list is printed with each number, as requested.
    print(*solvable_n_values, sep=", ")

find_solvable_hypercube_sizes()