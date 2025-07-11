import math

def find_min_N():
    """
    Finds the minimum positive integer N such that there exists a pair of
    positive integers (n1, n2) with n1 <= N and n2 <= N for which the group
    GG_{n1, n2}(r) is infinite.
    
    The condition for the group to be infinite is 1/n1 + 1/n2 <= 1.
    """
    N = 1
    while True:
        found_pair = None
        # Iterate through all pairs (n1, n2) where n1 <= N and n2 <= N.
        # To find the first solution efficiently, we can iterate from N downwards,
        # since larger n1, n2 make the sum 1/n1 + 1/n2 smaller.
        for n1 in range(1, N + 1):
            for n2 in range(1, N + 1):
                # The group is infinite if 1/n1 + 1/n2 <= 1.
                # Use a small tolerance for floating point comparison.
                if 1.0 / n1 + 1.0 / n2 <= 1.0 + 1e-9:
                    found_pair = (n1, n2)
                    break
            if found_pair:
                break
        
        if found_pair:
            n1, n2 = found_pair
            print(f"Searching for the smallest N...")
            print(f"For N = {N}, a pair (n1, n2) = ({n1}, {n2}) was found that makes the group infinite.")
            print("The condition for an infinite group is: 1/n1 + 1/n2 <= 1")
            
            # Create the final equation string
            # We explicitly output each number in the equation.
            sum_val = 1.0/n1 + 1.0/n2
            equation_str = f"1/{n1} + 1/{n2} = {sum_val}"

            print(f"Checking the condition for this pair: {equation_str}, which is indeed less than or equal to 1.")
            print(f"\nThus, the minimum positive integer N for which S(N) is not empty is {N}.")
            return N
            
        N += 1

if __name__ == "__main__":
    find_min_N()
