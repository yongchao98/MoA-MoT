import math

def calculate_norm(n):
    """
    Calculates the 1-norm of the correlation matrix T for a given even integer n.
    """
    if not isinstance(n, int) or n < 0 or n % 2 != 0:
        print("Error: Input 'n' must be a non-negative even integer.")
        return

    # The derived formula for the 1-norm is ||T||_1 = 2**(n+1) - (3**n - 1) / (3**n + 1)
    
    # Let's break down the final equation: ||T||_1 = A - B / C
    term1_A = 2**(n + 1)
    term2_numerator_B = 3**n - 1
    term2_denominator_C = 3**n + 1
    
    # The final result
    norm_value = term1_A - term2_numerator_B / term2_denominator_C

    print(f"For n = {n}, the 1-norm of the correlation matrix T is calculated using the formula:")
    print("||T||_1 = 2^(n+1) - (3^n - 1) / (3^n + 1)\n")
    
    print("Let's break down the calculation of each number in the equation:")
    print(f"The first term is: 2^({n}+1) = {term1_A}")
    print(f"The numerator of the fraction is: 3^{n} - 1 = {term2_numerator_B}")
    print(f"The denominator of the fraction is: 3^{n} + 1 = {term2_denominator_C}")
    print(f"The value of the fraction is {term2_numerator_B}/{term2_denominator_C} = {term2_numerator_B/term2_denominator_C}")
    
    print("\nPutting it all together, the final result is:")
    print(f"||T||_1 = {term1_A} - ({term2_numerator_B} / {term2_denominator_C}) = {norm_value}")

# --- Main execution ---
# You can change the value of n here to any non-negative even integer.
# For example, n=2, 4, 6, etc.
n = 4
calculate_norm(n)