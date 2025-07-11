def find_min_N():
    """
    Finds the minimum positive integer N for which there exists a pair
    of integers (n1, n2) such that n1 <= N, n2 <= N, and the group
    GG_{n1, n2}(r) is infinite for some r > 0.
    """
    # Start checking from N=2, as n=1 leads to finite groups.
    N = 2
    while True:
        # Flag to check if a valid pair is found for the current N
        found_pair = False
        solution_pair = None

        # Iterate through all possible pairs (n1, n2) for a given N
        # We start loops from 2 since n=1 gives a finite group.
        for n1 in range(2, N + 1):
            for n2 in range(2, N + 1):
                # The condition for the group to be infinite is 1/n1 + 1/n2 <= 1/2
                if 1 / n1 + 1 / n2 <= 0.5:
                    solution_pair = (n1, n2)
                    found_pair = True
                    break
            if found_pair:
                break
        
        # If a pair was found, this N is the minimum.
        if found_pair:
            print(f"Searching for the minimum positive integer N where S(N) is not empty.")
            print(f"The condition for an infinite group is 1/n₁ + 1/n₂ ≤ 1/2, with n₁, n₂ ≥ 2.")
            print(f"For N = {N-1}, no pair satisfied the condition.")
            print(f"For N = {N}, a valid pair (n₁, n₂) = {solution_pair} was found.")
            
            n1, n2 = solution_pair
            result = 1/n1 + 1/n2
            
            print(f"\nChecking the condition for this pair:")
            # As requested, output each number in the final equation.
            print(f"1/{n1} + 1/{n2} = {result}")
            print(f"Since {result} is less than or equal to 0.5, the condition is met.")
            print(f"Therefore, the minimum value for N is {N}.")
            break
        
        # If no pair was found, increment N and continue the search.
        N += 1

if __name__ == "__main__":
    find_min_N()
