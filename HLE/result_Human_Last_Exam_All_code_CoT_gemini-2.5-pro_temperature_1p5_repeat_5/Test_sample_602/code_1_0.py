import math

def calculate_l_n(n):
    """
    Calculates the exact value of l(n) using the derived analytical formula.
    
    The formula is derived from the step-by-step evaluation of the function:
    l(n) = f^(1)( f_M^(2)( f^(3)(P) ) )
    
    After simplification, the expression for l(n) is:
    l(n) = 2 - (2/n^2) * ((2n-1)*sqrt(n^2 - n + 1) - 1)
    """
    if not isinstance(n, int) or n < 5:
        raise ValueError("n must be an integer greater than or equal to 5.")

    # Apply the simplified analytical formula
    sqrt_term = math.sqrt(n**2 - n + 1)
    value = 2 - (2 / n**2) * ((2 * n - 1) * sqrt_term - 1)
    
    return value

def print_final_equation():
    """
    Prints the final symbolic equation for l(n).
    The numbers in the equation are 2, 2, 2, 1, 1, 2, 1, 1.
    """
    # Print the equation with its constituent numbers
    equation_string = "l(n) = 2 - (2/n**2) * ((2*n - 1)*sqrt(n**2 - n + 1) - 1)"
    print("The exact value of l(n) is given by the formula:")
    print(equation_string)

# Main execution
# As per the problem, n must be an integer >= 5.
# We will demonstrate the calculation for n=5.
n_val = 5
l_val = calculate_l_n(n_val)

# Print the formula first
print_final_equation()
print("\nFor a sample value of n = 5:")
print(f"l(5) = {l_val}")
# Return the calculated value for n=5 in the required format.
print(f"<<<{l_val}>>>")
