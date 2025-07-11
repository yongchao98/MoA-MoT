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
            r = n % 8
            if r == 3 or r == 5:
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
    if n == 1:
        return t
    else:
        return 0

def calculate_higher_central_charge(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for the U(1)_{2*N1} x U(1)_{-2*N2} theory.
    This formula is valid for odd n, N1, N2.
    """
    if n % 2 == 0:
        print(f"Warning: The formula for zeta_{n} is derived assuming n is odd.")
    if N1 % 2 == 0 or N2 % 2 == 0:
        print(f"Warning: The formula for zeta_{n} is derived assuming N1 and N2 are odd.")
    
    b = n
    c = N1 * N2
    
    # The Jacobi symbol (b/c) is what we need to compute.
    # We will use the law of quadratic reciprocity to compute (c/b) if needed,
    # but the derived formula is (n / N1*N2), which is (b/c).
    
    try:
        result = jacobi_symbol(b, c)
        print(f"For N1 = {N1}, N2 = {N2}, and n = {n}:")
        print(f"The higher central charge is given by the Jacobi symbol zeta_{n} = ({b} / {c})")
        
        # Intermediate steps for clarity as per the logic
        jacobi_N1 = jacobi_symbol(n, N1)
        jacobi_N2 = jacobi_symbol(n, N2)
        print(f"zeta_{n} = ({n} / {N1}) * ({n} / {N2}) = {jacobi_N1} * {jacobi_N2} = {jacobi_N1 * jacobi_N2}")
        print(f"Alternatively, zeta_{n} = ({n} / {N1 * N2}) = ({b} / {c}) = {result}")

    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    # Example values for N1, N2, and n.
    # Note: The formula is derived for odd N1, N2, and n.
    N1 = 3
    N2 = 5
    n = 7
    
    calculate_higher_central_charge(N1, N2, n)
    print("\n--- Another Example ---")
    N1 = 5
    N2 = 7
    n = 11
    calculate_higher_central_charge(N1, N2, n)
