import math

def calculate_phi(n):
    """
    Calculates the value of phi(n) based on the derived formula for n >= 5.

    The problem simplifies significantly under the given conditions. The final formula for phi(n) is
    exp(tr(Proj_M(X_inv))), where the trace of the projection is (n^2 - 1)*(2n - 3) / n.
    This script calculates and prints this value.
    """
    if not isinstance(n, int) or n < 5:
        print("Error: n must be an integer greater than or equal to 5.")
        return

    # The final expression for the exponent in phi(n) is (n**2 - 1) * (2*n - 3) / n.
    # As per the instructions, we output the components of this final equation.
    term1 = n**2 - 1
    term2 = 2 * n - 3
    denominator = n

    print(f"Calculating phi(n) for n = {n}")
    print("The final equation for the exponent is (n^2 - 1) * (2n - 3) / n")
    print(f"Value of '(n^2 - 1)': {term1}")
    print(f"Value of '(2n - 3)': {term2}")
    print(f"Value of 'n' (denominator): {denominator}")

    # Calculate the exponent (which is the trace of the projected matrix)
    exponent = (term1 * term2) / denominator
    print(f"\nThe exponent value is: {exponent}")

    # Calculate phi(n) = exp(exponent)
    result = math.exp(exponent)
    print(f"The final value of phi({n}) is: {result}")
    
    # Return the value for the final answer block
    return result

# The problem is defined for n >= 5. We will demonstrate the calculation for n=5.
n_value = 5
final_answer = calculate_phi(n_value)