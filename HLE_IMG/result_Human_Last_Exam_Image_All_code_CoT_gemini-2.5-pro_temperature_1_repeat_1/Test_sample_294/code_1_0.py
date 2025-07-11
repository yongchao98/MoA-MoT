import math

def find_min_N():
    """
    Finds the minimum positive integer N such that S(N) is non-empty.

    S(N) is the set of pairs (n1, n2) with n1, n2 <= N for which the
    group GG_{n1,n2}(r) is infinite for some r > 0.
    """
    print("The group GG_{n1,n2}(r) is infinite if and only if the corresponding")
    print("group of global plane rotations is infinite. This occurs if the condition")
    print("1/n1 + 1/n2 <= 1 is met. We also require n1 > 1 and n2 > 1, because if either")
    print("is 1, the group is trivially finite.")
    print("\nWe will search for the smallest integer N for which there exists a pair (n1, n2)")
    print("with 1 < n1 <= N and 1 < n2 <= N satisfying this inequality.\n")

    N = 1
    while True:
        # For a pair (n1, n2) to exist with n1, n2 > 1, N must be at least 2.
        # We can start the search loops from 2.
        if N >= 2:
            # Check all pairs (n1, n2) where 2 <= n1, n2 <= N
            for n1 in range(2, N + 1):
                for n2 in range(2, N + 1):
                    # Check if the condition for an infinite group is met.
                    # We use a small tolerance for floating point comparison.
                    if 1/n1 + 1/n2 <= 1.0000001:
                        print(f"Searching for N = {N}...")
                        print(f"Found the valid pair ({n1}, {n2}) where n1 <= {N} and n2 <= {N}.")
                        print("Checking the condition with this pair:")
                        
                        # Output the final equation with each number
                        result = 1/n1 + 1/n2
                        print(f"1 / {n1} + 1 / {n2} = {result}")

                        print(f"Since {result} <= 1, the condition for an infinite group is satisfied.")
                        print(f"Therefore, S({N}) is not empty.")
                        print(f"\nThe minimum value of N is {N}.")
                        return N
        N += 1

# Execute the function to find and print the result.
find_min_N()