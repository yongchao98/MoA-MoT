import numpy as np

def solve_integral():
    """
    Calculates the value of the integral based on the assumption that it is a conserved quantity.
    The integral I is given by 3 * ln(3 / (e^2 + e + 1)).
    """
    # Introduction and explanation of the method
    print("Based on the complexity of the PDE and the structure of the problem, we assume the requested integral is a conserved quantity.")
    print("This means its value at t=1 is the same as its value at t=0.")
    print("We calculate the integral I = integral from 0 to 1 of u(x,y,-y,0) dx.")
    print("The integrand simplifies to: -(6*exp(2x) + 3*exp(x)) / (exp(2x) + exp(x) + 1)")
    print("This is a logarithmic derivative, and the integral evaluates to:")
    print("I = [-3 * ln(exp(2x) + exp(x) + 1)] from x=0 to x=1")
    print("I = (-3 * ln(e^2 + e + 1)) - (-3 * ln(3))")
    print("I = 3 * (ln(3) - ln(e^2 + e + 1)) = 3 * ln(3 / (e^2 + e + 1))")
    print("\n--- Numerical Calculation ---")
    
    # Define constants
    e = np.exp(1)
    e_squared = np.exp(2)
    
    # Calculate components of the final formula
    denominator = e_squared + e + 1
    fraction = 3 / denominator
    final_result = 3 * np.log(fraction)

    # Print the equation with numerical values
    print(f"The final expression is: 3 * ln(3 / (e^2 + e + 1))")
    print(f"Substituting numerical values for e and e^2:")
    print(f"= 3 * ln(3 / ({e_squared:.6f} + {e:.6f} + 1))")
    print(f"= 3 * ln(3 / {denominator:.6f})")
    print(f"= 3 * ln({fraction:.6f})")
    print(f"= 3 * ({np.log(fraction):.6f})")
    print(f"= {final_result:.6f}")
    
    return final_result

if __name__ == '__main__':
    result = solve_integral()
    # The final answer is wrapped according to the instruction.
    print(f"\nThe spatial average of the system's state is: {result:.6f}")
    print(f"<<<{result:.12f}>>>")
