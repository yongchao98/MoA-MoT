import math

def solve_integral():
    """
    This function calculates the value of the definite integral based on the analytical solution.
    The integral is I = -3 * (ln(e^2 + e + 1) - ln(3)).
    """
    
    # Define the constants and terms in the final equation
    factor = -3
    e = math.e
    e_squared = e**2
    
    # Arguments for the natural logarithm
    log_arg1 = e_squared + e + 1
    log_arg2 = 3
    
    # Calculate the final result
    result = factor * (math.log(log_arg1) - math.log(log_arg2))
    
    # Print the equation and its components as requested
    print("The integral is calculated using the formula: I = C * (ln(A) - ln(B))")
    print("\nThe numbers and terms in this final equation are:")
    print(f"C (factor) = {factor}")
    print(f"e = {e}")
    print(f"A (e^2 + e + 1) = {e_squared} + {e} + 1 = {log_arg1}")
    print(f"B = {log_arg2}")
    
    # Print the step-by-step calculation
    print("\nCalculation steps:")
    print(f"ln(A) = ln({log_arg1}) = {math.log(log_arg1)}")
    print(f"ln(B) = ln({log_arg2}) = {math.log(log_arg2)}")
    print(f"I = {factor} * ({math.log(log_arg1)} - {math.log(log_arg2)})")
    print(f"I = {factor} * ({math.log(log_arg1) - math.log(log_arg2)})")
    
    # Print the final result
    print(f"\nThe final value of the integral is: {result}")

if __name__ == "__main__":
    solve_integral()
