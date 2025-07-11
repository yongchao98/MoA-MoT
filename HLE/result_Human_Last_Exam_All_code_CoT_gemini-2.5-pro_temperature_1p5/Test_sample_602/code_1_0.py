import math

def calculate_l_n(n):
    """
    Calculates the value of the function l(n) based on the derived analytical formula.
    
    The formula is: l(n) = (2*n^2 + 2 - 2*(2*n-1)*sqrt(n^2-n+1)) / n^2
    """
    if not isinstance(n, int) or n < 5:
        raise ValueError("n must be an integer greater than or equal to 5.")

    # Calculate the components of the formula
    # l(n) = (A - B * sqrt(C)) / D
    A = 2 * n**2 + 2
    B = 2 * (2 * n - 1)
    C = n**2 - n + 1
    D = n**2
    
    # Calculate the exact value
    sqrt_C = math.sqrt(C)
    numerator = A - B * sqrt_C
    result = numerator / D
    
    print(f"For n = {n}:")
    print(f"The final equation is of the form (A - B * sqrt(C)) / D")
    print(f"A = 2*n^2 + 2 = {A}")
    print(f"B = 2*(2*n - 1) = {B}")
    print(f"C = n^2 - n + 1 = {C}")
    print(f"D = n^2 = {D}")
    print(f"The exact symbolic value is ({A} - {B}*sqrt({C})) / {D}")
    print(f"The numerical value is: {result}")
    
    return result

# --- Main execution ---
# The user is asked to calculate the exact value of l(n) for n >= 5.
# As no specific n is given, we will use n=5 as a representative example.
n_value = 5 
final_value = calculate_l_n(n_value)

# The requested format is <<<value>>>.
# For n=5, the exact value is (52 - 18*sqrt(21))/25.
# The numerical value is approximately -1.21947.
# We will provide the numerical value as requested by the output format example.
print(f"\nFinal answer for n={n_value}:")
# print(f"<<<(52 - 18 * math.sqrt(21)) / 25>>>") # Exact Symbolic
# print(f"<<<{final_value}>>>") # Numerical
