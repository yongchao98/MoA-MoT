import math

def solve():
    """
    Calculates the minimum value of the sum of set sizes based on the given conditions.
    """
    try:
        n_str = input("Enter the value of n: ")
        n = int(n_str)
        if n <= 0:
            print("Please enter a positive integer for n.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # Based on the derivation, the sequence of set sizes |S_i| must be of the form x_i,
    # where x_i >= 1 and |x_{i+1} - x_i| = 1.
    # To minimize the sum, we should use the smallest possible integers,
    # which leads to an alternating sequence of 1s and 2s.

    # Sequence: 1, 2, 1, 2, 1, ...
    # If n is even, n = 2k, there are k ones and k twos. Sum = k * 1 + k * 2 = 3k = 3 * n / 2.
    # If n is odd, n = 2k + 1, there are k+1 ones and k twos. Sum = (k+1) * 1 + k * 2 = 3k + 1.
    # We can express this using integer division.
    
    num_of_pairs = n // 2
    sum_val = num_of_pairs * 3
    if n % 2 == 1:
        sum_val += 1
        
    # The derived values match choice B for n=2,3,4 but not for all n.
    # Choice B is floor(n^2/4) + 2. Let's calculate it as well.
    # Let's trust the logical derivation.

    # Final formula: 3 * (n // 2) + (n % 2)
    k = n // 2
    if n % 2 == 0:
        # For even n, the sum is 3*k. The sizes are 1, 2, 1, 2, ... , 1, 2
        min_sum = 3 * k
    else:
        # For odd n, the sum is 3*k + 1. The sizes are 1, 2, 1, 2, ..., 1
        min_sum = 3 * k + 1

    print(f"For n = {n}, the minimum value of the sum is: {min_sum}")

solve()