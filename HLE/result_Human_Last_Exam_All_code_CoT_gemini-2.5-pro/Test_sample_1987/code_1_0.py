import math
import sys

# Set a higher recursion limit for the recursive calculation, although it's
# not strictly necessary for the given small inputs like (2, 4, 5).
sys.setrecursionlimit(3000)

# Memoization table for the recursive function
memo = {}

def f_recursive(a):
    """
    Calculates the defined recursive function f(a_1, ..., a_n)
    using recursion with memoization to handle overlapping subproblems.
    """
    # Use a tuple for the state so it can be a dictionary key
    a = tuple(a)
    if a in memo:
        return memo[a]

    # Rule (2): Base case f(0, 0, ..., 0) = 1
    if all(x == 0 for x in a):
        return 1

    # Rule (1): Termination conditions
    # If a_1 < 0
    if a[0] < 0:
        return 0
    # If the sequence is not in increasing order
    for i in range(len(a) - 1):
        if a[i] > a[i+1]:
            return 0

    # Rule (3): Recursive step
    result = 0
    for i in range(len(a)):
        # Create the next state by subtracting 1 from the i-th element
        next_a = list(a)
        next_a[i] -= 1
        result += f_recursive(tuple(next_a))
    
    memo[a] = result
    return result

def solve():
    """
    This function calculates the answers to the three questions posed.
    """
    # --- Part 1: Calculate f(2, 4, 5) ---
    # The arguments (2, 4, 5) are small enough for the recursive approach.
    ans1 = f_recursive((2, 4, 5))

    # --- Part 2: Calculate f(9000, 9000, 9000) ---
    # For large arguments like (9000, 9000, 9000), recursion is infeasible.
    # We use the known closed-form formula for f(a, a, a):
    # f(a, a, a) = 2 * (3a)! / (a! * (a+1)! * (a+2)!)
    # To avoid overflow with huge factorials, we compute the logarithm of the value
    # using math.lgamma, where lgamma(x+1) = log(x!).
    a = 9000
    # log(f) = log(2) + log((3a)!) - log(a!) - log((a+1)!) - log((a+2)!)
    log_ans2 = (math.log(2) + math.lgamma(3 * a + 1) - 
                (math.lgamma(a + 1) + math.lgamma(a + 1 + 1) + math.lgamma(a + 2 + 1)))
    # Exponentiate to get the final value
    ans2 = math.exp(log_ans2)

    # --- Part 3: Calculate f(p, p, p, p) mod p for p = 10^9 + 7 ---
    # We use the formula for f(a,a,a,a) and evaluate it modulo p.
    # f(a,a,a,a) = C(4a;a,a,a,a) * (1!*2!*3!) / product_{k=1..3} [(a+1)...(a+k)]
    # Modulo p, the fractional part becomes 1.
    # The multinomial coefficient C(4p;p,p,p,p) mod p becomes C(4;1,1,1,1) by Lucas's Theorem.
    # C(4;1,1,1,1) = 4! / (1!*1!*1!*1!) = 24.
    p = 10**9 + 7
    ans3 = 24

    # Print the final answers, separated by a comma.
    print(f"{ans1},{ans2},{ans3}")

solve()