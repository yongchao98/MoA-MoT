import numpy as np
from scipy import integrate

def main():
    """
    Calculates the definite integral of a piecewise function from x=0 to x=4.
    """
    
    # Define the first part of the function for the interval [0, 3]
    def p1(x):
        return (2 * x**3) / 8

    # Define the second part of the function for the interval [3, 5]
    def p2(x):
        return np.exp(x) * (1 + np.sin(x)) / (1 + np.cos(x))

    # Calculate the integral of the first part from 0 to 3
    # integrate.quad returns a tuple (result, estimated_error)
    integral1, _ = integrate.quad(p1, 0, 3)

    # Calculate the integral of the second part from 3 to 4
    integral2, _ = integrate.quad(p2, 3, 4)

    # The total integral is the sum of the two parts
    total_integral = integral1 + integral2

    # Print the final equation with the calculated values
    print("The total integral is the sum of the integral from 0 to 3 and the integral from 3 to 4.")
    print(f"Integral part 1 (from 0 to 3): {integral1}")
    print(f"Integral part 2 (from 3 to 4): {integral2}")
    print("\nFinal Equation:")
    print(f"{integral1} + ({integral2}) = {total_integral}")

if __name__ == "__main__":
    main()