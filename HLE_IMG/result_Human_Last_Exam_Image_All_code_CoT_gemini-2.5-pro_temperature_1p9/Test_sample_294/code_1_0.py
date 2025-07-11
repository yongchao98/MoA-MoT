import math

def find_min_N():
    """
    Finds the minimum positive integer N for which S(N) is not empty.
    
    S(N) is the set of pairs (n1, n2) with n1, n2 <= N for which
    the group GG_{n1, n2}(r) is infinite for some r > 0.
    
    A known mathematical result states that for r=1, the group is infinite
    if and only if 1/n1 + 1/n2 <= 1/2. We search for the smallest N
    for which such a pair (n1, n2) exists.
    
    We only need to check n1, n2 >= 2, as for n=1 the corresponding
    rotation is the identity and the group is finite.
    """
    N = 2
    while True:
        found_pair = None
        # Iterate through all pairs (n1, n2) with 2 <= n1, n2 <= N
        for n1 in range(2, N + 1):
            for n2 in range(2, N + 1):
                # Check the condition for the group to be infinite
                if 1/n1 + 1/n2 <= 0.5:
                    found_pair = (n1, n2)
                    break
            if found_pair:
                break
        
        if found_pair:
            n1, n2 = found_pair
            print(f"Searching for N...")
            print(f"For N = {N}, we test the pair (n1, n2) = ({n1}, {n2}).")
            print("The condition for the group to be infinite is 1/n1 + 1/n2 <= 1/2.")
            print(f"Let's check for ({n1}, {n2}):")
            sum_of_reciprocals = 1/n1 + 1/n2
            print(f"1/{n1} + 1/{n2} = {1/n1:.2f} + {1/n2:.2f} = {sum_of_reciprocals:.2f}")
            if sum_of_reciprocals == 0.5:
              print(f"The sum is exactly 1/2, so the condition 1/n1 + 1/n2 <= 1/2 is met.")
            else:
              print(f"The sum is less than 1/2, so the condition 1/n1 + 1/n2 <= 1/2 is met.")

            print(f"\nSince a valid pair ({n1}, {n2}) exists for N = {N}, S({N}) is non-empty.")
            print(f"We have confirmed that for N < {N}, no such pair exists.")
            print(f"Therefore, the minimum N is {N}.")
            return N
        
        N += 1

find_min_N()