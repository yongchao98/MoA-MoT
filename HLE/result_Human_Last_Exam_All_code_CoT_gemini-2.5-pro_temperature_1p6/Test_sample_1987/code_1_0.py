import math
from functools import lru_cache
import sys

# Set a higher recursion limit for the first problem, although it might not be strictly necessary with good memoization
sys.setrecursionlimit(2000)

@lru_cache(maxsize=None)
def f(a: tuple) -> int:
    """
    Computes the recursive function f with memoization.
    The input 'a' must be a tuple to be hashable for the cache.
    """
    # Base case (1): a_1 < 0
    if a[0] < 0:
        return 0
    
    # Base case (1): sequence is not in increasing order
    for i in range(len(a) - 1):
        if a[i] > a[i+1]:
            return 0
            
    # Base case (2): f(0, 0, ..., 0) = 1
    if all(x == 0 for x in a):
        return 1
        
    # Recursive step (3)
    res = 0
    # Create a list from the tuple to modify its elements
    b_list = list(a)
    for i in range(len(a)):
        # Decrement the i-th element
        b_list[i] -= 1
        # Call the function recursively with the new tuple
        res += f(tuple(b_list))
        # Restore the original value for the next iteration
        b_list[i] += 1
    return res

def solve_part1():
    """
    Calculates f(2, 4, 5) using the recursive function.
    """
    return f((2, 4, 5))

def solve_part2():
    """
    Calculates f(9000, 9000, 9000) using its known combinatorial formula.
    f(k,k,k) = ( (3k)! / (k!)^3 ) * 2 / ((k+1)^2 * (k+2))
             = C(3k, k) * C(2k, k) * 2 / ((k+1)^2 * (k+2))
    """
    k = 9000
    
    # The formula can be computed using math.comb for binomial coefficients.
    # Python's integers support arbitrary size, so we can calculate this exactly.
    term1 = math.comb(3 * k, k)
    term2 = math.comb(2 * k, k)
    
    numerator = term1 * term2 * 2
    denominator = (k + 1)**2 * (k + 2)
    
    # The division is guaranteed to be exact.
    result = numerator // denominator
    return result

def solve_part3():
    """
    Calculates f(p, p, p, p) mod p for p = 10^9 + 7.
    The formula is f(k,k,k,k) = ( (4k)!/(k!)^4 ) * 6 / ( (k+1)^3 * (k+2)^2 * (k+3) )
    We evaluate this for k = p, where p is a prime number.
    """
    p = 10**9 + 7
    k = p
    
    # Numerator part ( (4k)!/(k!)^4 ) mod p:
    # This equals C(4k,k)*C(3k,k)*C(2k,k)*C(k,k) mod p.
    # By Lucas's Theorem C(ap,bp) = C(a,b) mod p.
    # So C(np,p) = C(n,1) = n mod p.
    # The term becomes 4 * 3 * 2 * 1 = 24.
    num_part1 = 24
    
    # Numerator part (6):
    num_part2 = 6
    
    # Denominator part ( (k+1)^3 * (k+2)^2 * (k+3) ) mod p:
    # For k=p, this is (p+1)^3 * (p+2)^2 * (p+3) mod p
    # which is congruent to 1^3 * 2^2 * 3 = 12 mod p.
    den_part = 12
    
    # The final result is (num_part1 * num_part2) / den_part mod p
    # Note: Division is modular inverse, but here we can simplify first.
    result = (num_part1 * num_part2) // den_part
    
    return result

def main():
    """
    Calculates and prints the results for the three parts of the problem.
    """
    ans1 = solve_part1()
    ans2 = solve_part2()
    ans3 = solve_part3()
    
    print(f"{ans1},{ans2},{ans3}")

main()