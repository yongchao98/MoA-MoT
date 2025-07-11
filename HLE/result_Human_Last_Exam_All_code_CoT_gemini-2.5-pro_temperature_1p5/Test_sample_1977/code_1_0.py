import math

def calculate_t_norm(n):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n.

    Args:
        n: A non-negative even integer.

    Returns:
        The value of the 1-norm. Returns None if n is invalid.
    """
    if not isinstance(n, int) or n < 0 or n % 2 != 0:
        print("Error: n must be a non-negative even integer.")
        return None

    # Numerator of the fraction part of the formula: 2**(n+2) + 4
    num = 4 * (2**n + 1)
    # Denominator of the fraction part of the formula: 3**n + 1
    den = 3**n + 1

    # The full formula: 2**(n+1) + 3 - num / den
    term1 = 2**(n+1) + 3
    result = term1 - num / den
    
    # We are asked to output each number in the final equation.
    # The equation is: result = (2**(n+1) + 3) - (4*(2**n+1))/(3**n+1)
    
    print(f"For n = {n}, the 1-norm of the correlation matrix T is calculated as:")
    print(f"||T||_1 = (2^({n}+1) + 3) - (4 * (2^{n} + 1)) / (3^{n} + 1)")
    
    # Calculate intermediate values for printing
    val_2_np1 = 2**(n+1)
    val_2_n_p1 = 2**n + 1
    val_4_mult = 4 * val_2_n_p1
    val_3_n_p1 = 3**n + 1
    term1_val = val_2_np1 + 3

    print(f"||T||_1 = ({val_2_np1} + 3) - (4 * {val_2_n_p1}) / ({val_3_n_p1})")
    print(f"||T||_1 = {term1_val} - {val_4_mult} / {val_3_n_p1}")

    # Display the final result
    if result.is_integer():
        print(f"Result: {int(result)}")
    else:
        # To display as a fraction
        common_divisor = math.gcd(val_4_mult, val_3_n_p1)
        num_simple = val_4_mult // common_divisor
        den_simple = val_3_n_p1 // common_divisor
        print(f"||T||_1 = {term1_val} - {num_simple}/{den_simple}")
        final_num = term1_val * den_simple - num_simple
        print(f"Result: {final_num}/{den_simple} (or approximately {result:.4f})")
    
    return result

if __name__ == '__main__':
    # You can change the value of n here.
    # Let's use n=2 as an example.
    n_value = 2
    final_answer = calculate_t_norm(n_value)
    
    # The problem asks for the formula for even n. The numerical result depends on the specific value of n.
    # Let's consider n=2 as a specific case.
    if n_value == 2 and final_answer is not None:
        print("\nFor n=2, the final answer is:")
        print(f"<<<{int(final_answer)}>>>")
