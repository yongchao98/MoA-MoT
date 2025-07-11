def find_min_N():
    """
    Finds the minimum positive integer N such that there exists a pair of
    positive integers (n1, n2) with n1 <= N and n2 <= N for which the
    group GG_{n1,n2}(r) is infinite.

    The condition for the group to be infinite is 1/n1 + 1/n2 <= 1.
    """
    N = 1
    while True:
        print(f"Checking for N = {N}...")
        found_pair = False
        # To avoid redundant checks, we can assume n1 <= n2
        for n1 in range(1, N + 1):
            for n2 in range(n1, N + 1):
                # We only need to check pairs that haven't been checked for smaller N
                if n1 < N and n2 < N:
                    continue

                value = 1.0 / n1 + 1.0 / n2
                print(f"  Testing pair ({n1}, {n2}): 1/{n1} + 1/{n2} = {value:.4f}")

                if value <= 1:
                    print(f"\nFound a valid pair ({n1}, {n2}) for N = {N}.")
                    print(f"The condition 1/n1 + 1/n2 <= 1 is met.")
                    # Final equation output as requested
                    print(f"Final equation check: 1/{n1} + 1/{n2} = {value}")
                    print(f"\nTherefore, the minimum N for which S(N) is not empty is {N}.")
                    return N
        
        print(f"S({N}) is empty.")
        N += 1

min_N = find_min_N()