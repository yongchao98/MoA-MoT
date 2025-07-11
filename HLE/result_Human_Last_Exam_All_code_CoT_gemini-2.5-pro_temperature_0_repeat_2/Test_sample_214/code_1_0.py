import collections
import sys

# It's good practice to increase recursion limit for deep poset structures.
sys.setrecursionlimit(2000)

def solve():
    """
    This function calculates the number of open sets in the topological space P^-(D_S, tau)
    by counting the number of order-preserving maps from the divisibility poset on S={1,...,150}
    to a 3-element chain.
    """
    N = 150
    
    # Find prime numbers up to N to help build the poset efficiently.
    primes = []
    is_prime = [True] * (N + 1)
    for p in range(2, N + 1):
        if is_prime[p]:
            primes.append(p)
            for i in range(p * p, N + 1, p):
                is_prime[i] = False

    # Build the poset by finding the children of each number.
    # A child of s is t=s*p where p is a prime and t <= N.
    children = collections.defaultdict(list)
    for s in range(1, N + 1):
        for p in primes:
            if s * p > N:
                break
            t = s * p
            children[s].append(t)
    
    maximal_elements = {s for s in range(1, N + 1) if not children[s]}

    # Memoization table for the dynamic programming.
    memo = {}

    def count_hom(s, c):
        """
        Recursively counts the number of valid colorings for the sub-poset
        of multiples of s, given that s has color c.
        """
        if (s, c) in memo:
            return memo[(s, c)]

        # Base case: For maximal elements, any color is valid.
        if s in maximal_elements:
            return 1

        # Recursive step: Product of sums of ways for children.
        total_ways = 1
        for child in children[s]:
            child_sum = 0
            # The child's color c_child must be >= parent's color c.
            for c_child in range(c, 3):
                child_sum += count_hom(child, c_child)
            total_ways *= child_sum
        
        memo[(s, c)] = total_ways
        return total_ways

    # The total number is the sum of possibilities for the root element '1'.
    n0 = count_hom(1, 0)
    n1 = count_hom(1, 1)
    n2 = count_hom(1, 2)
    total = n0 + n1 + n2

    print(f"The number of open sets is the sum of all valid order-preserving maps f: D_S -> C_3.")
    print(f"Number of maps with f(1)=0: {n0}")
    print(f"Number of maps with f(1)=1: {n1}")
    print(f"Number of maps with f(1)=2: {n2}")
    print(f"Total number of open sets = {n0} + {n1} + {n2} = {total}")
    
    return total

final_answer = solve()
print(f"<<<{final_answer}>>>")