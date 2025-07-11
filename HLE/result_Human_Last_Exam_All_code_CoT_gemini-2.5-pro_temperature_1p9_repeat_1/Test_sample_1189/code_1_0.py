import math
import sys

def solve_counting_problem():
    """
    Calculates the number of sets T satisfying the given conditions using a recurrence relation.

    The problem is to find the number of sets T such that:
    - |T| = m
    - T is a collection of non-empty subsets of S = {1, 2, ..., n}
    - For each i in S, i is contained in an even number of subsets in T.

    This is solved by translating the problem into finding the number of m-element
    subsets of non-zero vectors in the vector space (F_2)^n that sum to the zero vector.

    Let f(k) be this number for a set of size k. The function calculates f(m)
    using the recurrence relation:
    m * f(m) = C(N, m-1) - f(m-1) - (N - m + 2) * f(m-2)
    where N = 2^n - 1.

    The base cases are f(0) = 1 and f(1) = 0.
    """
    # Set the parameters n and m for the problem.
    # You can change these values to solve for different n and m.
    n = 4
    m = 4

    print(f"Solving for n = {n}, m = {m}\n")
    
    # The total number of non-empty subsets of S
    N = 2**n - 1
    
    # If m is larger than the number of available subsets, the answer is 0.
    if m > N:
        print("Since m is greater than the total number of non-empty subsets (2^n - 1),")
        print("it's impossible to choose m distinct subsets.")
        print("f({m}) = 0")
        print("\n<<<0>>>")
        return

    # A list to store the results of f(k) for k = 0, 1, ..., m
    f = [0] * (m + 1)
    
    # Base case f(0) = 1 (the empty set's sum is zero)
    if m >= 0:
        f[0] = 1
        print("Base case: f(0) = 1")
    
    # Base case f(1) = 0 (a single non-zero vector can't sum to zero)
    if m >= 1:
        f[1] = 0
        print("Base case: f(1) = 0")
    
    # Iteratively compute f(k) for k from 2 to m
    for k in range(2, m + 1):
        # Term 1: Combinations C(N, k-1)
        # math.comb(n, k) returns 0 if k > n
        term1 = math.comb(N, k - 1)
        
        # Term 2: f(k-1)
        term2 = f[k - 1]
        
        # Term 3: (N - k + 2) * f(k-2)
        factor3 = N - k + 2
        term3 = factor3 * f[k - 2]
        
        # Numerator of the recurrence
        numerator = term1 - term2 - term3
        
        # Result for f(k)
        result = numerator // k
        f[k] = result
        
        # Print the calculation steps for the current k
        print(f"\nCalculating f({k}):")
        print(f"{k} * f({k}) = C({N}, {k-1}) - f({k-1}) - ({N} - {k} + 2) * f({k-2})")
        print(f"{k} * f({k}) = {term1} - {term2} - {factor3} * {f[k-2]}")
        print(f"{k} * f({k}) = {term1} - {term2} - {term3}")
        print(f"{k} * f({k}) = {numerator}")
        print(f"f({k}) = {result}")

    final_answer = f[m]
    print(f"\nFinal Answer: The number of such sets T is {final_answer}.")
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    solve_counting_problem()