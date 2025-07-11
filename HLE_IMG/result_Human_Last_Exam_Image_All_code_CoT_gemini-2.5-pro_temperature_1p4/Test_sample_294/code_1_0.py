def find_min_N():
    """
    Finds the minimum positive integer N for which S(N) is non-empty.
    
    The condition for the group GG_{n1, n2}(r) to be infinite for some r > 0
    is n1 >= 2 and n2 >= 2.
    
    The code iterates through N = 1, 2, 3, ... and for each N, it checks if
    there exists a pair (n1, n2) such that 1 <= n1, n2 <= N and the condition
    for infinitude is met.
    """
    N = 1
    while True:
        s_is_non_empty = False
        found_pair = None
        # Check if S(N) is non-empty for the current N
        for n1 in range(1, N + 1):
            for n2 in range(1, N + 1):
                # Condition for the group to be infinite
                if n1 >= 2 and n2 >= 2:
                    s_is_non_empty = True
                    found_pair = (n1, n2)
                    break
            if s_is_non_empty:
                break
        
        if s_is_non_empty:
            print(f"For N = {N}, a pair (n1, n2) satisfying the conditions for an infinite group exists.")
            print(f"The first such pair is n1 = {found_pair[0]}, n2 = {found_pair[1]}.")
            print("This pair satisfies the required inequalities:")
            # Showing the "equation" part as requested
            print(f"n1 = {found_pair[0]} >= 2")
            print(f"n2 = {found_pair[1]} >= 2")
            print(f"\nSince N = {N} is the smallest positive integer for which such a pair exists,")
            print(f"the minimum value for N is {N}.")
            return N
        
        N += 1

# Execute the search and print the final result.
min_N_value = find_min_N()