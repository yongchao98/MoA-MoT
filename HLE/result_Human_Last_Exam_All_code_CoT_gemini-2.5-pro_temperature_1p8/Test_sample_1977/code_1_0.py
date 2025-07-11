import math

def calculate_norm_for_J_n(n):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n with an even n.

    Args:
        n (int): An even non-negative integer.

    Returns:
        float: The 1-norm of the correlation matrix.
    """
    if not isinstance(n, int) or n < 0 or n % 2 != 0:
        raise ValueError("n must be a non-negative even integer.")

    term1_val = 2**(n + 1)
    term2_val = 3
    
    num_val = 4 * (2**n + 1)
    den_val = 3**n + 1
    
    term3_val = num_val / den_val
    
    result = term1_val + term2_val - term3_val
    
    print(f"For n = {n}:")
    print(f"The formula is: 2**(n+1) + 3 - (4 * (2**n + 1)) / (3**n + 1)")
    print(f"Substituting n = {n}:")
    print(f"2**({n+1}) + 3 - (4 * (2**{n} + 1)) / (3**{n} + 1)")
    print(f"= {term1_val} + {term2_val} - (4 * ({2**n} + 1)) / ({3**n} + 1)")
    print(f"= {term1_val} + {term2_val} - ({num_val}) / ({den_val})")
    print(f"= {term1_val + term2_val} - {term3_val}")
    print(f"= {result}")

if __name__ == '__main__':
    # You can change the value of n here. It must be a non-negative even integer.
    # For example, let's calculate for n=2 and n=4.
    try:
        print("--- Calculation for n=2 ---")
        calculate_norm_for_J_n(2)
        print("\n" + "="*20 + "\n")
        print("--- Calculation for n=4 ---")
        calculate_norm_for_J_n(4)
    except ValueError as e:
        print(e)
