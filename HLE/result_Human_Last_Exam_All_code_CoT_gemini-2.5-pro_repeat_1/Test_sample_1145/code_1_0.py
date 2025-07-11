import numpy as np
from scipy.integrate import quad

def solve_probability():
    """
    Calculates the probability of finding a particle in a 1D box
    within a specified interval for a given quantum state n.
    """
    # --- Problem Parameters ---
    # Quantum state
    n = 2
    # Length of the box. We can set a=1.0 since the interval is a
    # fraction of 'a', making the result independent of 'a'.
    a = 1.0
    # Interval boundaries
    x1 = 0.495 * a
    x2 = 0.505 * a

    # --- Print the Equation Being Solved ---
    # The probability P is the integral of the probability density |psi(x)|^2
    # over the interval [x1, x2].
    # P = Integral from x1 to x2 of (2/a) * sin^2(n*pi*x/a) dx
    print("This script calculates the probability by solving the following integral:")
    print(f"P = Integral from x = {x1:.3f}a to x = {x2:.3f}a of (2/a) * sin^2({n} * pi * x / a) dx")
    print("-" * 60)

    # --- Define the Integrand ---
    # This is the probability density function for the particle.
    def probability_density(x, n_val, a_val):
        """
        |psi_n(x)|^2 for a particle in a 1D box.
        """
        return (2.0 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

    # --- Perform Numerical Integration ---
    # The quad function returns the result of the integral and an error estimate.
    # We pass the constant parameters n and a using the 'args' argument.
    probability, error = quad(probability_density, x1, x2, args=(n, a))

    # --- Display the Final Result ---
    print(f"The calculated probability for n = {n} in the interval [{x1}, {x2}] is:")
    # We use scientific notation for the small result.
    print(f"{probability:.5e}")

    # Return the value for the final answer block
    return probability

# Run the calculation and print the results
final_answer = solve_probability()
# The final answer tag requires a direct value.
# We format the calculated probability to be placed there.
# For example, <<<0.000006580>>> or <<<6.580e-6>>>.
# Let's use the decimal representation.
final_answer_str = f"<<<{final_answer:.10f}>>>"
# This final print is for the thought process and will not be in the final block.
# To conform to the format, I will manually place the calculated value at the end.

if __name__ == '__main__':
    # This block is not strictly necessary for the user, but it's good practice.
    # The function call above already prints the output.
    pass
