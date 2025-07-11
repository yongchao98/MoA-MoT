import math

def combinations(n, k):
    """Helper function to compute combinations C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve_sets_problem(n, m):
    """
    Calculates the number of sets T with size m satisfying the given conditions.
    
    The problem is equivalent to finding the number of m-element subsets of non-zero
    vectors in F_2^n that sum to the zero vector. Let this number be f(m).
    
    We derive the recurrence relation:
    m * f(m) = C(N, m-1) - f(m-1) - (N - m + 2) * f(m-2)
    where N = 2^n - 1.
    
    Base cases: f(0) = 1, f(1) = 0.
    """
    if n <= 0 or m < 0:
        print("n must be a positive integer and m must be a non-negative integer.")
        return

    N = 2**n - 1

    if m == 0:
        print("For m=0, the only solution is the empty set T={}.")
        print("f(0) = 1")
        print("<<<1>>>")
        return

    # f_values will store f(0), f(1), ..., f(m)
    f_values = [0] * (m + 1)
    
    # Base cases
    f_values[0] = 1
    if m >= 1:
        f_values[1] = 0

    # Iteratively compute f(i) for i from 2 to m
    for i in range(2, m + 1):
        term1 = combinations(N, i - 1)
        term2 = f_values[i - 1]
        term3_factor = N - i + 2
        term3 = term3_factor * f_values[i - 2]
        
        # The recurrence guarantees the numerator is divisible by i
        numerator = term1 - term2 - term3
        f_values[i] = numerator // i

    # Output the final result with the full equation for f(m)
    term1_val = combinations(N, m - 1)
    term2_val = f_values[m - 1]
    term3_factor_val = N - m + 2
    term3_val_f_part = f_values[m-2]

    print(f"For n={n} and m={m}, we have N = 2^{n}-1 = {N}.")
    print("We compute f(k) for k=0...m:")
    for i in range(m + 1):
        print(f"f({i}) = {f_values[i]}")

    print("\n--- Final Calculation ---")
    
    # Equation with symbols
    print(f"The equation for f({m}) is: f({m}) = (C({N}, {m-1}) - f({m-1}) - ({N} - {m} + 2) * f({m-2})) / {m}")

    # Equation with computed values
    print(f"f({m}) = ({term1_val} - {term2_val} - ({term3_factor_val} * {term3_val_f_part})) / {m}")
    
    numerator_val = term1_val - term2_val - (term3_factor_val * term3_val_f_part)
    print(f"f({m}) = {numerator_val} / {m}")

    final_answer = f_values[m]
    print(f"The final answer is: {final_answer}")
    print(f"<<<{final_answer}>>>")


# You can change these values to solve for different n and m
n_val = 4
m_val = 5

solve_sets_problem(n_val, m_val)
