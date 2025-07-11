import math

def jacobi(a, n):
    """
    Computes the Jacobi symbol (a/n).
    Assumes n is a positive odd integer.
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

def calculate_zeta_n(n, N1, N2):
    """
    Calculates the higher central charge zeta_n for U(1)_{2N1} x U(1)_{-2N2} theory.
    """
    print(f"Calculating zeta_{n} for N1 = {N1}, N2 = {N2}")
    
    # Check if N1 and N2 are positive odd integers
    if N1 <= 0 or N1 % 2 == 0:
        print("Error: N1 must be a positive odd integer.")
        return
    if N2 <= 0 or N2 % 2 == 0:
        print("Error: N2 must be a positive odd integer.")
        return

    # Components of the formula: zeta_n = (2/N1)*(n/N1)*(-1/N2)*(2/N2)*(n/N2)
    # The term (2n/N1)
    j_2n_N1 = jacobi(2 * n, N1)
    
    # The term (-2n/N2)
    j_neg2n_N2 = jacobi(-2 * n, N2)
    
    # The final result is the product of these two
    zeta_n_val = j_2n_N1 * j_neg2n_N2

    # Breaking down the terms for detailed output
    j_2_N1 = jacobi(2, N1)
    j_n_N1 = jacobi(n, N1)
    j_neg1_N2 = jacobi(-1, N2)
    j_2_N2 = jacobi(2, N2)
    j_n_N2 = jacobi(n, N2)
    
    print("\nThe formula is: zeta_n = (2n/N1) * (-2n/N2)")
    print("Which expands to: (2/N1) * (n/N1) * (-1/N2) * (2/N2) * (n/N2)")

    print("\nIntermediate values of the Jacobi symbols:")
    print(f"(2n/N1) = ({2*n}/{N1}) = {j_2n_N1}")
    print(f"  (2/N1) = (2/{N1}) = {j_2_N1}")
    print(f"  (n/N1) = ({n}/{N1}) = {j_n_N1}")

    print(f"(-2n/N2) = ({-2*n}/{N2}) = {j_neg2n_N2}")
    print(f"  (-1/N2) = (-1/{N2}) = {j_neg1_N2}")
    print(f"  (2/N2) = (2/{N2}) = {j_2_N2}")
    print(f"  (n/N2) = ({n}/{N2}) = {j_n_N2}")

    print(f"\nFinal equation with numbers:")
    print(f"zeta_{n} = ({j_2n_N1}) * ({j_neg2n_N2}) = {zeta_n_val}")
    print(f"  = ({j_2_N1})*({j_n_N1}) * ({j_neg1_N2})*({j_2_N2})*({j_n_N2})")
    
    # Returning the final answer to be wrapped
    return zeta_n_val


if __name__ == '__main__':
    # Example usage:
    # Please enter integer values for n, N1, and N2.
    # Note: N1 and N2 must be positive odd integers for the Jacobi symbol to be well-defined.
    try:
        n_in = int(input("Enter integer n: "))
        N1_in = int(input("Enter positive odd integer N1: "))
        N2_in = int(input("Enter positive odd integer N2: "))
        
        result = calculate_zeta_n(n_in, N1_in, N2_in)
        if result is not None:
            # The final answer as requested
            print(f"\n<<< {result} >>>")

    except ValueError:
        print("\nInvalid input. Please enter integers.")
