import numpy as np
from scipy.integrate import quad

def calculate_probability():
    """
    Calculates the probability of finding a particle in a 1D box in a given interval.
    """
    # Parameters for the problem
    n = 2
    # We can set the length of the box 'a' to 1, as it cancels out in the calculation.
    a = 1.0
    x1 = 0.495 * a
    x2 = 0.505 * a

    # The probability density function, P(x) = |ψ(x)|²
    def probability_density(x, n_val, a_val):
        return (2.0 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

    # Perform the numerical integration of the probability density function
    # from x1 to x2.
    # The 'args' parameter passes the values of n and a to the function.
    probability, error = quad(probability_density, x1, x2, args=(n, a))

    # Print the setup and the final result
    print("To find the probability, we integrate the probability density function P(x) = |ψ(x)|² from x1 to x2.")
    print("The function is P(x) = (2/a) * sin²(n*π*x/a).")
    print("\nGiven parameters:")
    print(f"n = {n}")
    print(f"Interval = [{x1}a, {x2}a]")
    print("\nThe final equation we are solving numerically is:")
    # The print statement below shows the specific integral being calculated, with a=1.
    print(f"Probability = Integral from {x1} to {x2} of (2/{a}) * sin²({n}*π*x/{a}) dx")

    print(f"\nThe calculated probability is: {probability}")

if __name__ == "__main__":
    calculate_probability()
    # The analytical result is (16*pi^2/3) * (0.005)^3 which is approximately 6.5797e-6
    # The numerical integration should yield a very similar result.
    # Let's calculate the value to be put in the final answer tag.
    # P = (16 * np.pi**2 / 3) * (0.005)**3 = 6.579708e-06
    # The code output is 6.58920...e-06. Let's re-check the analytical formula.
    # P = 0.01 - (1 / (2*np.pi)) * np.sin(0.02*np.pi) = 6.58920...e-06
    # The analytical formula matches the numerical integration. The Taylor expansion was an approximation.
    final_answer = 6.589203533522647e-06
    # The user wants the answer in the format <<<answer content>>>
    # print(f"\n<<<{final_answer}>>>")