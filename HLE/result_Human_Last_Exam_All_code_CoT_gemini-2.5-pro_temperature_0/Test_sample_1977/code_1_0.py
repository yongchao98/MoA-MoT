import math

def calculate_t_norm_1(n: int):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n.

    Args:
        n: An even integer representing the parameter of the state J_n.
    """
    if n % 2 != 0 or n < 0:
        print("Error: n must be a non-negative even integer.")
        return

    print(f"Calculating the 1-norm of the correlation matrix T for n = {n}.")
    print("-" * 30)

    # The derived formula for the 1-norm is:
    # ||T||_1 = 2^(n+1) + 3 - 4 * (2^n + 1) / (3^n + 1)
    print("The formula for the 1-norm is:")
    print("||T||_1 = 2**(n+1) + 3 - 4 * (2**n + 1) / (3**n + 1)")
    print("\nSubstituting n =", n)

    # Calculate term by term
    term1_val = 2**(n + 1) + 3
    print(f"\nStep 1: Calculate the first part of the expression.")
    print(f"2**(n+1) + 3 = 2**({n+1}) + 3 = {2**(n+1)} + 3 = {term1_val}")

    numerator = 4 * (2**n + 1)
    denominator = 3**n + 1
    term2_val = numerator / denominator
    
    print(f"\nStep 2: Calculate the second part (the fraction).")
    print(f"4 * (2**n + 1) / (3**n + 1) = 4 * (2**{n} + 1) / (3**{n} + 1)")
    print(f"= 4 * ({2**n} + 1) / ({3**n} + 1)")
    print(f"= 4 * {2**n + 1} / {3**n + 1}")
    print(f"= {numerator} / {denominator}")
    
    # Check if the result is an integer
    if numerator % denominator == 0:
        print(f"= {int(term2_val)}")
    else:
        # Simplify fraction if possible
        common_divisor = math.gcd(numerator, denominator)
        simplified_num = numerator // common_divisor
        simplified_den = denominator // common_divisor
        if common_divisor > 1:
            print(f"= {simplified_num}/{simplified_den} ≈ {term2_val:.4f}")
        else:
            print(f"≈ {term2_val:.4f}")


    # Final result
    result = term1_val - term2_val
    
    print("\nStep 3: Combine the parts to get the final result.")
    print(f"||T||_1 = {term1_val} - {term2_val}")
    
    if isinstance(result, int) or result.is_integer():
        print(f"||T||_1 = {int(result)}")
    else:
        print(f"||T||_1 = {result}")
        
    print("-" * 30)
    return result

# --- Main execution ---
# You can change the value of n here to any non-negative even integer.
# For example, n = 0, 2, 4, 6, ...
n = 2
final_answer = calculate_t_norm_1(n)

# Example for another even n
# n = 4
# calculate_t_norm_1(n)
