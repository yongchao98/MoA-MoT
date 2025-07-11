import math

def find_min_N():
    """
    Finds the minimum positive integer N for which S(N) is non-empty.
    """
    # We are looking for the minimum N such that there exists a pair (n1, n2)
    # with 1 <= n1, n2 <= N for which the group GG_{n1,n2}(r) is infinite.

    # The condition for the group to be infinite is n1 >= 2, n2 >= 2, and
    # cos(pi/n1) + cos(pi/n2) < 1.

    # We can iterate through N starting from 1 and find the first N for which
    # such a pair (n1, n2) exists.
    N = 1
    while True:
        found_pair = False
        # For a given N, search for a valid pair (n1, n2)
        for n1 in range(1, N + 1):
            for n2 in range(1, N + 1):
                # Condition 1: n1 and n2 must be at least 2.
                if n1 < 2 or n2 < 2:
                    continue

                # Condition 2: cos(pi/n1) + cos(pi/n2) < 1
                value = math.cos(math.pi / n1) + math.cos(math.pi / n2)

                if value < 1:
                    # A valid pair is found for this N. Since we are iterating N
                    # upwards, this N is the minimum.
                    print(f"Searching for the minimum N...")
                    print(f"For N = {N}, we check for pairs (n1, n2) with 2 <= n1, n2 <= {N}.")
                    print(f"The first pair that can be tested is ({n1}, {n2}) when N reaches {N}.")
                    print("\nThe condition for the group to be infinite is: cos(pi/n1) + cos(pi/n2) < 1")
                    print(f"Testing the pair (n1, n2) = ({n1}, {n2}):")
                    
                    cos_n1 = math.cos(math.pi / n1)
                    cos_n2 = math.cos(math.pi / n2)
                    
                    print(f"cos(pi/{n1}) + cos(pi/{n2}) = {cos_n1:.4f} + {cos_n2:.4f} = {value:.4f}")
                    print(f"Since {value:.4f} is less than 1, the condition is satisfied.")
                    print(f"Therefore, S({N}) is not empty.")
                    print(f"\nAs N=1 does not allow any pairs with n1,n2 >= 2, the minimum value for N is {N}.")
                    return N
        
        N += 1

min_N_value = find_min_N()
print(f"<<<{min_N_value}>>>")
