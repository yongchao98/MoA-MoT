import sys
import math

# The recursion depth can be large for the first part, so we increase the limit.
# For f(2,4,5), the max depth is 2+4+5=11, which is fine, but this is a good practice.
sys.setrecursionlimit(2000)

# Memoization dictionary for the first calculation (f(2,4,5))
memo_f1 = {}

def f_recursive(a_tuple):
    """
    This function calculates f(a_1, a_2, ..., a_n) based on the provided recursive definition.
    It uses memoization (dynamic programming) to store and retrieve results for previously 
    computed states, which is crucial for performance.
    """
    # Base case: check if result is already memoized
    if a_tuple in memo_f1:
        return memo_f1[a_tuple]

    # Rule (1): f(...) = 0 if a_1 < 0.
    # The recursion never generates a_i < 0 for i > 1 if the list is sorted.
    if a_tuple[0] < 0:
        return 0
    
    # Rule (1): f(...) = 0 if the sequence is not in increasing order
    is_sorted = all(a_tuple[i] <= a_tuple[i+1] for i in range(len(a_tuple) - 1))
    if not is_sorted:
        return 0
        
    # Rule (2): f(0, 0, ..., 0) = 1
    if all(x == 0 for x in a_tuple):
        return 1

    # Rule (3): The recursive step
    # f(a) = sum_i f(a with a_i decremented)
    res = 0
    for i in range(len(a_tuple)):
        next_a_list = list(a_tuple)
        next_a_list[i] -= 1
        res += f_recursive(tuple(next_a_list))
        
    # Memoize the result before returning
    memo_f1[a_tuple] = res
    return res

# Part 1: Calculate f(2, 4, 5)
# We use the recursive function with memoization, which works well for small inputs.
val1_tuple = (2, 4, 5)
val1 = f_recursive(val1_tuple)

# Part 2: Calculate f(9000, 9000, 9000)
# A direct recursive calculation would be too slow. The value f(k, k, ..., k) (n times) 
# is given by the Fuss-Catalan formula, A_k^(n) = (1/((n-1)k+1)) * C(nk, k).
# For n=3 and k=9000, this is f(9000,9000,9000) = (1/(2*9000+1)) * C(3*9000, 9000).
k_val2 = 9000
n_val2 = 3
# We use Python's math.comb for combinations, which handles large integers.
# The Fuss-Catalan numbers are integers, so the division is exact.
val2 = math.comb(n_val2 * k_val2, k_val2) // ((n_val2 - 1) * k_val2 + 1)

# Part 3: Calculate f(p, p, p, p) mod p for p=10^9+7
# We use the same Fuss-Catalan formula with n=4 and k=p.
# We need to compute A_p^(4) mod p = (1/((4-1)p+1)) * C(4p, p) mod p.
# p is the prime number 10^9+7.
# Denominator: (3p+1) === 1 (mod p). Its modular inverse is 1.
# Numerator: C(4p, p) mod p is evaluated using Lucas's Theorem.
# Base-p expansion of 4p is (4, 0) and p is (1, 0).
# So, C(4p, p) === C(4, 1) * C(0, 0) (mod p).
# C(4, 1) = 4 and C(0, 0) = 1.
# Thus, C(4p, p) === 4 (mod p).
# The final result is (1)^-1 * 4 === 4 (mod p).
p = 10**9 + 7
val3 = 4

print(f"{val1},{val2},{val3}")