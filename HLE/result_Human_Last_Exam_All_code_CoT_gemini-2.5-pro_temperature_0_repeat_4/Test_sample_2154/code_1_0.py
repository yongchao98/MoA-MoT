import math

def gcd(a, b):
    """
    Computes the greatest common divisor of two integers using the Euclidean algorithm.
    """
    while b:
        a, b = b, a % b
    return a

def count_new_fractions(k):
    """
    Calculates the number of irreducible fractions j/k such that 0 < j/k < 1/2.
    This corresponds to the number of integers j such that 1 <= j < k/2 and gcd(j, k) = 1.
    """
    if k <= 2:
        return 0
    
    count = 0
    # We need j/k < 1/2, which means j < k/2.
    # The loop for j goes from 1 up to, but not including, k/2.
    # math.ceil(k / 2) gives the correct upper bound for the range.
    limit_j = math.ceil(k / 2)
    for j in range(1, limit_j):
        if gcd(j, k) == 1:
            count += 1
    return count

def solve_sequence():
    """
    Calculates and prints the sequence u_r(n) for n from 3 to 12.
    """
    results = {}
    
    print("Calculating the sequence u_r(n) for n = 3 to 12:")
    
    # Base case: For n=3, k can be 1 or 2. There are no fractions satisfying the conditions.
    ur_val = 0
    results[3] = ur_val
    print(f"u_r(3) = {ur_val}")

    # Iteratively calculate u_r(n) for n = 4 to 12
    for n in range(4, 13):
        # The value for u_r(n) is u_r(n-1) plus the contribution from k = n-1.
        k = n - 1
        new_contrib = count_new_fractions(k)
        
        prev_ur_val = ur_val
        ur_val = prev_ur_val + new_contrib
        results[n] = ur_val
        
        print(f"u_r({n}) = u_r({n-1}) + (contribution from k={k}) = {prev_ur_val} + {new_contrib} = {ur_val}")

    print("\n--------------------------------------------------")
    print("The complete set of values is:")
    result_list = [results[n] for n in range(3, 13)]
    print(f"{{u_r(3), u_r(4), ..., u_r(12)}} = {result_list}")

# Execute the solution
solve_sequence()