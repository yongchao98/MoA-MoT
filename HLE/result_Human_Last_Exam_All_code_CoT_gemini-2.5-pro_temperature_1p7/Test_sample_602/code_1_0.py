import numpy as np

def calculate_l(n: int):
    """
    Calculates the exact value of l(n) using the derived formula.

    The formula for l(n) is:
    l(n) = 2 + 2/n^2 - (4n-2)/n^2 * sqrt(n^2-n+1)

    This script will print the steps to evaluate the formula for the given n
    and output the final numerical result.
    """
    if not isinstance(n, int) or n < 5:
        raise ValueError("Input n must be an integer greater than or equal to 5.")

    # Calculate components of the formula
    n_squared = n**2
    term1 = 2
    term2_num = 2
    term2_den = n_squared
    
    term3_num_coeff = 4 * n - 2
    term3_den = n_squared
    
    sqrt_term_val = n**2 - n + 1
    sqrt_val = np.sqrt(sqrt_term_val)

    # Construct the equation string
    print(f"To calculate l(n) for n = {n}, we use the derived formula:")
    print(f"l(n) = 2 + 2/n^2 - (4n-2)/n^2 * sqrt(n^2-n+1)")
    print("\nSubstituting n = {}:".format(n))
    print(f"l({n}) = 2 + 2/({n}^2) - (4*{n}-2)/({n}^2) * sqrt({n}^2-{n}+1)")
    print(f"l({n}) = {term1} + {term2_num}/{term2_den} - {term3_num_coeff}/{term3_den} * sqrt({sqrt_term_val})")

    # Final calculation
    result = term1 + term2_num / term2_den - (term3_num_coeff / term3_den) * sqrt_val
    
    print("\nFinal Result:")
    print(f"l({n}) = {result}")
    
    return result

if __name__ == '__main__':
    # As per the problem description, n >= 5. Let's calculate for n=5.
    n_value = 5
    final_value = calculate_l(n_value)
    # The final answer format is specified to be <<<value>>>
    # print(f"\n<<<{final_value}>>>")