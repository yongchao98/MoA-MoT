import math

def find_min_N():
    """
    This function finds the minimum positive integer N for which S(N) is not empty.
    S(N) is the set of pairs (n1, n2) with n1, n2 <= N for which the group GG_{n1,n2}(r) is infinite.
    The condition for the group to be infinite is 1/n1 + 1/n2 <= 1, and we must have n1, n2 >= 2.
    """
    N = 1
    while True:
        s_n_is_non_empty = False
        found_n1 = -1
        found_n2 = -1
        
        # We need to find if there exists a pair (n1, n2) for the current N
        # We search from N downwards to find the largest integers first, which are more likely to satisfy the condition
        for n1 in range(N, 0, -1):
            for n2 in range(N, 0, -1):
                # The group is finite if n1=1 or n2=1. So we require n1, n2 >= 2.
                # The condition for the group to be infinite is 1/n1 + 1/n2 <= 1.
                if n1 >= 2 and n2 >= 2:
                    if (1/n1 + 1/n2) <= 1:
                        s_n_is_non_empty = True
                        found_n1 = n1
                        found_n2 = n2
                        break
            if s_n_is_non_empty:
                break
        
        if s_n_is_non_empty:
            print(f"Searching for the minimal N...")
            print(f"For N = {N}, a valid pair (n1, n2) = ({found_n1}, {found_n2}) was found.")
            print(f"This pair requires N to be at least {max(found_n1, found_n2)}.")
            print("\nThe condition for the corresponding group to be infinite is satisfied:")
            
            result = 1/found_n1 + 1/found_n2
            
            # Print the equation with all the numbers
            print(f"1 / {found_n1} + 1 / {found_n2} = {result}")
            print(f"Since {result} <= 1, the condition holds.")
            
            print(f"\nThus, the minimum positive integer N for which S(N) is non-empty is {N}.")
            print(f"<<<{N}>>>")
            return
            
        N += 1

# Execute the function to find the answer.
find_min_N()
