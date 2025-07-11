import math

def find_min_N():
    """
    Finds the minimum positive integer N such that S(N) is not empty.

    S(N) is the set of pairs (n1, n2) with n1, n2 <= N for which the group
    GG_{n1, n2}(r) is infinite for some r > 0.
    
    The condition for the group to be infinite is 1/n1 + 1/n2 <= 1,
    with n1, n2 >= 2.
    """
    print("Starting the search for the minimum N...")
    
    N = 1
    while True:
        print(f"\n--- Checking N = {N} ---")
        found_pair = None
        
        # Iterate through all possible pairs (n1, n2) where n1 <= N and n2 <= N
        for n1 in range(1, N + 1):
            for n2 in range(1, N + 1):
                # The order of rotations must be at least 2 for non-trivial rotation.
                # If n1=1 or n2=1, the group is finite.
                if n1 < 2 or n2 < 2:
                    continue

                # Check if the pair (n1, n2) satisfies the condition for an infinite group.
                if (1 / n1) + (1 / n2) <= 1:
                    found_pair = (n1, n2)
                    break
            if found_pair:
                break
        
        if found_pair:
            n1, n2 = found_pair
            print(f"Found a valid pair (n1, n2) = ({n1}, {n2}) for N = {N}.")
            print(f"This pair satisfies the conditions: n1 <= {N}, n2 <= {N}.")
            print("The condition for the group to be infinite is 1/n1 + 1/n2 <= 1.")
            
            val1 = 1 / n1
            val2 = 1 / n2
            result = val1 + val2
            
            print("\nVerification of the condition for the pair (", n1, ", ", n2, "):")
            print(f"1/{n1} + 1/{n2} = {val1} + {val2} = {result}")
            
            if result <= 1:
                print(f"Since {result} <= 1, the condition is met.")
                print(f"\nThis means the set S(N) for N={N} is not empty.")
                print(f"Therefore, the minimum value for N is {N}.")
                return N
        else:
            print(f"No pair (n1, n2) found for N = {N} that satisfies the conditions.")
            print("S(N) is empty. Trying the next N.")
            N += 1

if __name__ == '__main__':
    min_N = find_min_N()
    # The final answer is the integer value found.
