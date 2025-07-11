import math

def solve(n, m):
    """
    Calculates the number of sets T of size m satisfying the given conditions.

    The problem is equivalent to finding the number of m-element subsets of
    non-zero vectors in the vector space F_2^n that sum to the zero vector.

    Let f_m be this number. The solution is based on the following recurrence relation:
    f_0 = 1
    For m > 0:
    - If m is odd:  f_m = (1/m) * C(2^n - 1, m - 1) - f_{m-1}
    - If m is even: f_m = ((2^n - m)/m) * f_{m-1} which simplifies to
                     f_m = ((2^(n-1) - m/2) / (m/2)) * f_{m-1} for integer arithmetic.

    Args:
        n (int): A positive integer for the set S = {1, 2, ..., n}.
        m (int): A positive integer for the size of the set T.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n <= 0 or m < 0:
        print("Error: n must be a positive integer and m must be a non-negative integer.")
        return

    if m == 0:
        print("f_0 = 1")
        print("The number of such sets is 1 (the empty collection of subsets).")
        print("\n<<<1>>>")
        return

    # f_i stores the value for f_i
    # Initialize with f_0 = 1
    f_prev = 1

    # Iteratively compute f_i from i=1 to m
    for i in range(1, m + 1):
        if i % 2 == 1:  # i is odd
            comb_val = math.comb(2**n - 1, i - 1)
            # f_i = comb_val / i - f_{i-1}
            # The formula m*(f_m + f_{m-1}) = comb(...) ensures the division is exact.
            f_curr = comb_val // i - f_prev
            if i == m:
                print(f"f_{i} = C(2^{n}-1, {i-1}) / {i} - f_{i-1}")
                print(f"f_{i} = C({2**n-1}, {i-1}) / {i} - {f_prev}")
                print(f"f_{i} = {comb_val} / {i} - {f_prev}")
                print(f"f_{i} = {comb_val // i} - {f_prev} = {f_curr}")

        else:  # i is even
            k = i // 2
            # f_i = ((2**(n-1) - k) / k) * f_{i-1}
            # The formula k*f_{2k} = (2**(n-1)-k)*f_{2k-1} ensures the division is exact.
            numerator = (2**(n - 1) - k) * f_prev
            f_curr = numerator // k
            if i == m:
                print(f"f_{i} = (2^({n}-1) - {i}//2) * f_{i-1} / ({i}//2)")
                print(f"f_{i} = ({2**(n-1)} - {k}) * {f_prev} / {k}")
                print(f"f_{i} = {2**(n-1) - k} * {f_prev} / {k}")
                print(f"f_{i} = {numerator} / {k} = {f_curr}")
        
        f_prev = f_curr

    final_answer = f_prev
    print(f"\nThe number of such sets T for n={n}, m={m} is {final_answer}.")
    print(f"\n<<<{final_answer}>>>")

# Example usage:
# You can change the values of n and m to solve for different inputs.
if __name__ == '__main__':
    # Example from the problem description context
    n_val = 4
    m_val = 5
    solve(n_val, m_val)