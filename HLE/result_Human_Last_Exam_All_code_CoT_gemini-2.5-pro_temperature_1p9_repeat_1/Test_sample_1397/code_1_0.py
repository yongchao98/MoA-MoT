import math

def is_prime(num):
    """Checks if a number is prime."""
    if num < 2:
        return False
    for i in range(2, int(math.sqrt(num)) + 1):
        if num % i == 0:
            return False
    return True

def is_composite(num):
    """Checks if a number is composite."""
    if num < 4:
        return False
    return not is_prime(num)

def find_graph_parameters(n_max):
    """
    Attempts to find integer solutions for the number of vertices belonging to
    0, 1, or 2 C5-cycles based on the problem's constraints.

    The constraints lead to the following system of equations:
    1) n_0 + n_1 + n_2 = n
    2) n_1 + 2*n_2 = 5n
    where n_0, n_1, n_2 are non-negative integers.

    Algebraic manipulation of these equations yields:
    n_2 = 4*n + n_0
    
    This function will demonstrate that no non-negative integer solutions
    for n_0, n_1, n_2 exist for any valid n.
    """
    print("Searching for the smallest composite n such that a graph with the given properties exists...")
    print("The properties imply the following relations for n, n_0, n_1, n_2:")
    print("1) n_0 + n_1 + n_2 = n")
    print("2) n_1 + 2*n_2 = 5n")
    print("where n_k is the number of vertices belonging to k C5s.\n")

    found_solution = False
    # A 7-regular graph must have n >= 8 and n must be even.
    for n in range(8, n_max + 1):
        if not is_composite(n) or n % 2 != 0:
            continue

        print(f"--- Checking for n = {n} ---")
        solution_for_n = False
        # Based on algebra (n_2 = 4*n + n_0), we know n_2 must be > n.
        # But for completeness, we can iterate through all possibilities to show none work.
        for n_0 in range(n + 1):
            for n_2 in range(n + 1):
                # From n_2 - n_0 = 4n
                if n_2 - n_0 == 4 * n:
                    # We have a candidate solution for n_0 and n_2.
                    # Let's find n_1 from n_0 + n_1 + n_2 = n
                    n_1 = n - n_0 - n_2
                    if n_1 >= 0:
                        # Now we must check if it satisfies the second original equation
                        # n_1 + 2*n_2 = 5*n
                        print(f"Checking potential solution: n_0={n_0}, n_1={n_1}, n_2={n_2}")
                        print(f"Equation 1: {n_0} + {n_1} + {n_2} = {n_0+n_1+n_2} (must be {n})")
                        print(f"Equation 2: {n_1} + 2*({n_2}) = {n_1 + 2*n_2} (must be {5*n})")

                        if n_0 + n_1 + n_2 == n and n_1 + 2*n_2 == 5*n:
                            print(f"Found a valid solution for n = {n}: (n_0, n_1, n_2) = ({n_0}, {n_1}, {n_2})")
                            solution_for_n = True
                            found_solution = True
                            break
            if solution_for_n:
                break
        
        if not solution_for_n:
            print(f"No solution found for n = {n}. As derived from the equations, n_2 must be >= 4*n, but n_2 cannot be larger than n.")
        
        if found_solution:
            break
        print("")

    if not found_solution:
        print("\nSearch complete.")
        print("Conclusion: The analytical derivation is correct. No such graph exists for any n > 0 because the problem's conditions are self-contradictory.")

if __name__ == '__main__':
    find_graph_parameters(50)
