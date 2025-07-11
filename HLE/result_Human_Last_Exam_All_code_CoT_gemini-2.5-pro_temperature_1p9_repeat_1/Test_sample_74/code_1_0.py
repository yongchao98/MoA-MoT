import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
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
    if n == 1:
        return t
    else:
        return 0

def calculate_higher_central_charge(n, N1, N2):
    """
    Calculates the higher central charge zeta_n for U(1)_{2*N1} x U(1)_{-2*N2}.
    This formula is valid for odd n.
    """
    if n % 2 == 0:
        return "Undefined for even n"
        
    try:
        # Jacobi symbol ( (N1*N2)/n )
        jacobi_N1N2_n = jacobi_symbol(N1 * N2, n)
        
        # Jacobi symbol ( -1/n ) which is (-1)^((n-1)/2)
        jacobi_neg1_n = 1 if (n - 1) // 2 % 2 == 0 else -1

        zeta_n = jacobi_neg1_n * jacobi_N1N2_n
        
        print(f"For n={n}, N1={N1}, N2={N2}:")
        print(f"The Jacobi symbol (-1 / {n}) is: {jacobi_neg1_n}")
        print(f"The Jacobi symbol (({N1}*{N2}) / {n}) is: {jacobi_N1N2_n}")
        print(f"zeta_{n} = (-1/{n}) * (N1*N2/{n}) = {jacobi_neg1_n} * {jacobi_N1N2_n} = {zeta_n}")
        return zeta_n
        
    except ValueError as e:
        return str(e)


# Example values from a physical context.
# n must be an odd integer for the formula and Jacobi symbol to be well-defined.
n_val = 3
N1_val = 5
N2_val = 7

result = calculate_higher_central_charge(n_val, N1_val, N2_val)
print(f"\nFinal calculated value for zeta_{n_val} is: {result}")
<<<1>>>