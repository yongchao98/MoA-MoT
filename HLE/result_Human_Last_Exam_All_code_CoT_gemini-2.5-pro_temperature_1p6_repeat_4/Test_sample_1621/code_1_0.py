import math

def get_t(n):
    """Calculates the exponent t in n = s * 2^t, where s is odd."""
    if n == 0:
        return 0
    t = 0
    while n % 2 == 0:
        n //= 2
        t += 1
    return t

def hurwitz_radon_rho(n):
    """
    Calculates the Hurwitz-Radon function rho(n).
    Let n = s * 2^t where s is odd.
    Let t = 4a + b, where 0 <= b <= 3.
    Then rho(n) = 8a + 2^b.
    """
    if n <= 0:
        return 0
    t = get_t(n)
    a = t // 4
    b = t % 4
    return 8 * a + (1 << b) # 2**b

def solve_problem():
    """
    Finds the natural numbers n for which the condition holds.
    Based on a known mathematical result, this is true if and only if
    n belongs to the set {1, 2, 4, 8}.
    We can find these values by checking the inequality n <= rho(n),
    where rho(n) is the Hurwitz-Radon function.
    """
    possible_n = []
    # We check n up to a reasonable limit to demonstrate the pattern.
    # The mathematical proof confirms there are no other solutions.
    limit = 100 
    for n in range(1, limit + 1):
        if n <= hurwitz_radon_rho(n):
            possible_n.append(n)
            
    print(f"The set of natural numbers n for which such matrices exist is: {possible_n}")
    print(f"Thus, there are {len(possible_n)} such natural numbers.")

# Run the solver
solve_problem()
