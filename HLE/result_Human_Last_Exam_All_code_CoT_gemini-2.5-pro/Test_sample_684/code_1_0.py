def find_possible_n_values():
    """
    This function finds all integers n in the range [7, 55] for which it is possible
    to reach a state with a single gift on an n x ... x n hypercube.
    The condition for this is that n must be congruent to 1 or 13 modulo 14.
    """
    
    possible_n = []
    print("Checking values of n from 7 to 55...\n")
    
    for n in range(7, 56):
        # The problem is solvable if and only if n = 14k + 1 or n = 14k + 13
        # for some non-negative integer k.
        
        if n % 14 == 1:
            k = (n - 1) // 14
            print(f"For n = {n}: It is a solution, since {n} = 14 * {k} + 1.")
            possible_n.append(n)
        elif n % 14 == 13:
            k = (n - 13) // 14
            print(f"For n = {n}: It is a solution, since {n} = 14 * {k} + 13.")
            possible_n.append(n)
            
    print("\nThe possible values for n in increasing order are:")
    # The problem asks to output the numbers, so we print the final list.
    print(', '.join(map(str, possible_n)))

# Execute the function
find_possible_n_values()
