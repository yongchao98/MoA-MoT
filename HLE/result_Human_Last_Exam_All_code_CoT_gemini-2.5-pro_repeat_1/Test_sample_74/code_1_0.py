import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("The denominator n must be a positive odd integer.")
    if a == 0:
        return 1 if n == 1 else 0
    if a == 1:
        return 1
    
    a = a % n
    if a == 0:
        return 1 if n == 1 else 0

    # Handle factors of 2 in a
    t = 1
    while a % 2 == 0:
        a //= 2
        n_mod_8 = n % 8
        if n_mod_8 == 3 or n_mod_8 == 5:
            t = -t
    
    # Quadratic Reciprocity
    if a % 4 == 3 and n % 4 == 3:
        t = -t
    
    # Now a is odd, so we can flip
    return t * jacobi_symbol(n, a)

def calculate_zeta_n(n, N1, N2):
    """
    Calculates the higher central charge zeta_n for U(1)_{2N1} x U(1)_{-2N2}.
    """
    print(f"Calculating zeta_{n} for n={n}, N1={N1}, N2={N2}")

    if n % 2 == 0:
        print("The formula is derived for odd n. The result for even n may not be well-defined.")
        return

    # Step 1: Calculate gcds
    d1 = math.gcd(n, N1)
    d2 = math.gcd(n, N2)
    print(f"d1 = gcd({n}, {N1}) = {d1}")
    print(f"d2 = gcd({n}, {N2}) = {d2}")

    # Step 2: Calculate arguments for Jacobi symbols
    n_d1 = n // d1
    N1_d1 = N1 // d1
    n_d2 = n // d2
    N2_d2 = N2 // d2

    # Check if denominators for Jacobi symbol are valid
    if N1_d1 <= 0 or N1_d1 % 2 == 0:
        print(f"Error: Denominator N1/d1 = {N1_d1} is not a positive odd integer. Jacobi symbol is not defined.")
        j1 = "undefined"
        result = "undefined"
    else:
        j1 = jacobi_symbol(n_d1, N1_d1)

    if N2_d2 <= 0 or N2_d2 % 2 == 0:
        print(f"Error: Denominator N2/d2 = {N2_d2} is not a positive odd integer. Jacobi symbol is not defined.")
        j2 = "undefined"
        result = "undefined"
    else:
        j2 = jacobi_symbol(n_d2, N2_d2)

    # Step 3: Compute the product
    if j1 != "undefined" and j2 != "undefined":
        result = j1 * j2
    else:
        result = "undefined"

    # Step 4: Print the final equation
    print(f"The formula for zeta_{n} is:")
    print(f"zeta_{n} = ( (n/d1) / (N1/d1) ) * ( (n/d2) / (N2/d2) )")
    print(f"zeta_{n} = ( {n_d1} / {N1_d1} ) * ( {n_d2} / {N2_d2} )")
    print(f"zeta_{n} = {j1} * {j2}")
    print(f"Final Result: zeta_{n} = {result}")
    return result

if __name__ == '__main__':
    # Example values, can be changed by the user.
    n = 15
    N1 = 35
    N2 = 21

    final_answer = calculate_zeta_n(n, N1, N2)
    # The final answer in the required format
    if final_answer != "undefined":
        print(f"\n<<< {final_answer} >>>")
