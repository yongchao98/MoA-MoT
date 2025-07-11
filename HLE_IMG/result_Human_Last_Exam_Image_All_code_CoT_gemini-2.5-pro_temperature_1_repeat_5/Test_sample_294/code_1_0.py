import math

def find_min_N():
    """
    This function finds the minimum positive integer N such that the set S(N) is non-empty.

    S(N) is the set of pairs of positive integers (n1, n2) such that:
    1. n1 <= N and n2 <= N
    2. There exists some r > 0 for which the group GG_{n1,n2}(r) is infinite.

    The condition for the group to be infinite for some r > 1 is:
    1/2 < 1/n1 + 1/n2 < 1
    
    The rotations are trivial for n=1, so we only consider n1, n2 >= 2.
    """
    N = 1
    while True:
        # We start searching for pairs from n1, n2 = 2
        # because for n=1, the rotation is trivial and the group is finite.
        if N < 2:
            N += 1
            continue

        found_pair = False
        # Iterate through all pairs (n1, n2) such that n1 <= N and n2 <= N
        for n1 in range(2, N + 1):
            for n2 in range(2, N + 1):
                # The sum of reciprocals
                s = 1/n1 + 1/n2
                
                # Check the condition for the group to be infinite
                if 0.5 < s < 1:
                    print(f"Searching for N...")
                    print(f"Found a valid pair for N = {N}: (n1, n2) = ({n1}, {n2})")
                    print("\nThis pair satisfies the condition for the group to be infinite for some r > 1.")
                    print("The condition is: 1/2 < 1/n1 + 1/n2 < 1")
                    print(f"For our pair ({n1}, {n2}), the equation is:")
                    # Output each number in the final equation
                    print(f"1/2 < 1/{n1} + 1/{n2} < 1")
                    print(f"0.5 < {s:.4f} < 1")
                    print("This inequality is true.")
                    
                    print(f"\nSince a valid pair was found for N={N}, and no pair exists for N < {N},")
                    print(f"the minimum integer N is {N}.")
                    found_pair = True
                    break
            if found_pair:
                break
        
        if found_pair:
            return N
        
        N += 1

# Run the function to find and print the answer
min_N = find_min_N()
print(f"\nFinal Answer: {min_N}")