import math

def calculate_correlation_norm(n: int):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n with an even n.

    The formula for the 1-norm is derived as:
    ||T||_1 = (2**(n+1) * 3**n + 3**(n+1) - 2**(n+1) - 1) / (1 + 3**n)
    
    Args:
        n: A non-negative even integer.
    """
    if n < 0 or n % 2 != 0:
        print("Error: n must be a non-negative even integer.")
        return

    # The formula for the 1-norm ||T||_1 is S / Z, where:
    # Z is the normalization factor of the state J_n.
    # S is the sum of the numerators of all |t_mu| terms.
    
    # 1. Calculate the normalization factor Z = 1 + 3^n
    Z_val = 1 + 3**n

    # 2. Calculate the total sum S based on the derived formula:
    # S = 2^(n+1) * 3^n + 3^(n+1) - 2^(n+1) - 1
    term1 = 2**(n+1) * 3**n
    term2 = 3**(n+1)
    term3 = -(2**(n+1))
    term4 = -1
    S_val = term1 + term2 + term3 + term4

    # 3. Calculate the final 1-norm
    norm_T = S_val / Z_val

    # 4. Output the calculation step-by-step
    print(f"Calculating the 1-norm of the correlation matrix for n = {n}:")
    print("-" * 50)
    
    print("The final formula for the 1-norm is ||T||_1 = S / Z, where:")
    print("Z = 1 + 3^n")
    print("S = 2^(n+1) * 3^n + 3^(n+1) - 2^(n+1) - 1")
    print("\nStep-by-step calculation:")
    
    # Equation for Z
    print(f"1. Z = 1 + 3^{n} = 1 + {3**n} = {Z_val}")

    # Equation for S
    print("2. S = 2^(n+1) * 3^n + 3^(n+1) - 2^(n+1) - 1")
    print(f"   S = 2^({n+1}) * 3^{n} + 3^({n+1}) - 2^({n+1}) - 1")
    # Print numbers in the equation for S
    print(f"   S = {2**(n+1)} * {3**n} + {3**(n+1)} - {2**(n+1)} - 1")
    # Print evaluated terms
    print(f"   S = {term1} + {term2} + ({term3}) + ({term4})")
    print(f"   S = {S_val}")

    # Final result
    print("\n3. ||T||_1 = S / Z")
    print(f"   ||T||_1 = {S_val} / {Z_val}")
    print(f"   ||T||_1 = {norm_T}")

    # To fulfill the final answer format request
    # print(f"\n<<<{norm_T}>>>")

if __name__ == '__main__':
    # You can change the value of n here to any non-negative even integer.
    # For example, n = 0, 2, 4, 6, ...
    even_n = 4
    calculate_correlation_norm(even_n)
