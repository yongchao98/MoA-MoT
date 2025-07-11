import math

def find_min_N():
    """
    This function finds the minimum positive integer N for which there exists
    a pair of positive integers (n1, n2) such that n1 <= N, n2 <= N,
    and the condition 1/n1 + 1/n2 < 1 is met.
    """
    print("Searching for min{N} such that S(N) is not empty.")
    print("The condition for a pair (n1, n2) to be in S(N) is 1/n1 + 1/n2 < 1.")
    print("-" * 40)
    
    N = 1
    while True:
        print(f"Testing N = {N}...")
        found_pair = False
        
        # Iterate through all pairs (n1, n2) where 1 <= n1, n2 <= N
        for n1 in range(1, N + 1):
            for n2 in range(1, N + 1):
                # Check if the pair satisfies the "hyperbolic" condition
                if 1 / n1 + 1 / n2 < 1:
                    print(f"Found a valid pair (n1, n2) = ({n1}, {n2}) for N = {N}.")
                    # Output the equation with the numbers, as requested.
                    sum_val = 1/n1 + 1/n2
                    print(f"Condition check: 1/{n1} + 1/{n2} = {sum_val:.4f}, which is < 1.")
                    print("-" * 40)
                    print(f"The smallest integer N for which S(N) is non-empty is {N}.")
                    found_pair = True
                    # Once a pair is found, we have our minimal N and can stop.
                    return N
        
        if not found_pair:
            print(f"For N = {N}, no pair satisfied the condition.")
            N += 1

# Execute the function and print the final answer in the required format.
min_N = find_min_N()
print(f"<<<{min_N}>>>")
