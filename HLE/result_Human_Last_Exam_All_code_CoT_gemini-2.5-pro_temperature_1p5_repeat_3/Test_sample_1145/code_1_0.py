import numpy as np

def calculate_probability():
    """
    Calculates the probability of finding a particle in a 1D box in a given
    interval for a specific quantum number n.
    """
    # Define the quantum state and the interval boundaries (normalized by box length 'a')
    n = 2
    u1 = 0.495  # Start of the interval (x1/a)
    u2 = 0.505  # End of the interval (x2/a)

    # The probability is calculated from the analytical solution of the integral of |Ψ(x)|².
    # The formula is: P = (u2 - u1) - [sin(2*n*pi*u2) - sin(2*n*pi*u1)] / (2*n*pi)

    # Calculate the terms of the formula
    interval_width = u2 - u1
    
    arg1 = 2 * n * np.pi * u1
    arg2 = 2 * n * np.pi * u2
    
    sin_term = np.sin(arg2) - np.sin(arg1)
    
    denominator = 2 * n * np.pi

    # Calculate the final probability
    probability = interval_width - (sin_term / denominator)

    # As requested, print the final equation with each number explicitly shown.
    # The 'pi' symbol is used for clarity instead of its numerical value.
    print("The final equation with all numbers plugged in is:")
    print(f"P = ({u2} - {u1}) - [sin(2 * {n} * pi * {u2}) - sin(2 * {n} * pi * {u1})] / (2 * {n} * pi)")

    # Print the final numerical result
    print("\nThe calculated probability is:")
    print(probability)

if __name__ == '__main__':
    calculate_probability()
<<<6.579785860475345e-06>>>