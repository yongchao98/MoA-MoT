import math
import cmath

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    
    a %= n
    if a == 0:
        return 0
    if a == 1:
        return 1
    
    # Rule (2/n) = (-1)^((n^2-1)/8)
    if a == 2:
        if n % 8 == 1 or n % 8 == 7:
            return 1
        else: # n % 8 == 3 or n % 8 == 5
            return -1

    # Law of Quadratic Reciprocity
    if a >= n:
        return jacobi_symbol(a % n, n)
    elif a % 2 == 0:
        return jacobi_symbol(2, n) * jacobi_symbol(a // 2, n)
    else: # a is odd and smaller than n
        if a % 4 == 3 and n % 4 == 3:
            return -jacobi_symbol(n, a)
        else:
            return jacobi_symbol(n, a)

def calculate_higher_central_charge(n, N1, N2):
    """
    Calculates the higher central charge zeta_n for the theory U(1)_{2*N1} x U(1)_{-2*N2}.
    This formula is valid for odd integers n.
    """
    print(f"Calculating for n={n}, N1={N1}, N2={N2}")
    
    if n % 2 == 0:
        print("Error: This formula is valid only for odd n.")
        return

    # Step 1: Calculate d1 = gcd(n, N1) and d2 = gcd(n, N2)
    d1 = math.gcd(n, N1)
    d2 = math.gcd(n, N2)
    print(f"d1 = gcd({n}, {N1}) = {d1}")
    print(f"d2 = gcd({n}, {N2}) = {d2}")

    # The final formula for zeta_n is:
    # zeta_n = exp(i*pi*n*(d1-d2)/(4*d1*d2)) * ( (N1/d1) / (n/d1) ) * ( (N2/d2) / (n/d2) )
    # where (a/b) is the Jacobi symbol.

    # Step 2: Calculate the exponential part
    exp_numerator = 1j * math.pi * n * (d1 - d2)
    exp_denominator = 4 * d1 * d2
    exp_part = cmath.exp(exp_numerator / exp_denominator)
    
    print("\n--- Final Equation Components ---")
    print("zeta_n = (Exponential Part) * (Jacobi Symbol 1) * (Jacobi Symbol 2)")
    print(f"Exponential Part: exp(i*pi*{n}*({d1}-{d2}) / (4*{d1}*{d2})) = {exp_part}")

    # Step 3: Calculate the first Jacobi symbol
    n_prime1 = n // d1
    N_prime1 = N1 // d1
    if math.gcd(N_prime1, n_prime1) != 1:
         # This case should not happen based on the definition of gcd
         print(f"Warning: gcd(N1/d1, n/d1) = gcd({N_prime1}, {n_prime1}) is not 1.")

    jacobi1 = jacobi_symbol(N_prime1, n_prime1)
    print(f"Jacobi Symbol 1: ({N_prime1}/{n_prime1}) = {jacobi1}")

    # Step 4: Calculate the second Jacobi symbol
    n_prime2 = n // d2
    N_prime2 = N2 // d2
    if math.gcd(N_prime2, n_prime2) != 1:
        # This case should not happen
        print(f"Warning: gcd(N2/d2, n/d2) = gcd({N_prime2}, {n_prime2}) is not 1.")
        
    jacobi2 = jacobi_symbol(N_prime2, n_prime2)
    print(f"Jacobi Symbol 2: ({N_prime2}/{n_prime2}) = {jacobi2}")
    
    # Step 5: Calculate the final result
    zeta_n = exp_part * jacobi1 * jacobi2
    
    print("\n--- Final Result ---")
    print(f"zeta_{n} = {exp_part} * {jacobi1} * {jacobi2} = {zeta_n}")
    
    # For comparison, print the symbolic formula
    formula = f"zeta_{n} = exp(i*pi*n*(d1-d2)/(4*d1*d2)) * ( (N1/d1) / (n/d1) ) * ( (N2/d2) / (n/d2) )"
    print("\nSymbolic Formula:")
    print(formula)
    print("where d1=gcd(n,N1), d2=gcd(n,N2), and (a/b) is the Jacobi symbol.")


if __name__ == '__main__':
    # Example 1: n=3, N1=5, N2=7
    # In this case d1=1, d2=1, so the exponential term is 1.
    # zeta_3 = (5/3)*(7/3) = (-1)*(1) = -1
    calculate_higher_central_charge(3, 5, 7)
    
    print("\n" + "="*40 + "\n")
    
    # Example 2: n=3, N1=6, N2=5
    # d1=3, d2=1. Exp part is exp(i*pi*3*(3-1)/(4*3*1)) = exp(i*pi/2) = i
    # Jacobi1 = (6/3 / 3/3) = (2/1) = 1
    # Jacobi2 = (5/1 / 3/1) = (5/3) = -1
    # zeta_3 = i * 1 * (-1) = -i
    calculate_higher_central_charge(3, 6, 5)

<<<For an odd integer $n$, the higher central charge $\zeta_n$ for the theory $U(1)_{2N_1} \times U(1)_{-2N_2}$ is given by the formula:
$$
\zeta_n = e^{\frac{i\pi n(d_1-d_2)}{4d_1d_2}} \left(\frac{N_1/d_1}{n/d_1}\right) \left(\frac{N_2/d_2}{n/d_2}\right)
$$
where $d_1 = \gcd(n, N_1)$, $d_2 = \gcd(n, N_2)$, and $\left(\frac{a}{b}\right)$ is the Jacobi symbol.>>>