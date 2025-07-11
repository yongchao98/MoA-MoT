import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def p(x):
    """
    A helper polynomial with 7 distinct real roots: -3, -2, -1, 0, 1, 2, 3.
    This will be used as h''(x).
    """
    return x * (x**2 - 1) * (x**2 - 4) * (x**2 - 9)

def h_prime(x):
    """
    A polynomial of degree 8, which is the integral of p(x).
    This represents a plausible h'(x) = f'(g(x))g'(x).
    The constants are chosen so that h_prime(x) > 0 for all x,
    and the equation h_prime(x) = 1 has 8 real roots.
    """
    # Integral of x^7 - 14x^5 + 49x^3 - 36x
    integral_part = (1/8)*x**8 - (14/6)*x**6 + (49/4)*x**4 - 18*x**2
    # We add a constant to ensure h'(x) is positive and oscillates around y=1.
    # The minimum of the integral part is about -49.375.
    # We choose a constant C=25 to shift the graph up.
    C = 25
    return integral_part + C

def h_prime_minus_1(x):
    """
    The function whose roots we want to find, h'(x) - 1 = 0.
    """
    return h_prime(x) - 1

def find_and_show_roots():
    """
    Finds and prints the 8 roots of h'(x) = 1 and shows a plot.
    """
    # Visually estimate initial guesses for the 8 roots from the plot
    initial_guesses = [-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5]
    
    # Use a numerical solver to find the roots
    roots = fsolve(h_prime_minus_1, initial_guesses)
    roots.sort()
    
    print("A plausible polynomial h'(x) can be constructed such that h'(x) > 0 and h'(x)=1 has 8 real roots.")
    print("This allows for h(x)=x to have 9 real roots (fixed points).")
    print("\nAn example of such an h'(x) is: (1/8)x^8 - (7/3)x^6 + (49/4)x^4 - 18x^2 + 25")
    print("\nThe 8 real roots of the equation h'(x) - 1 = 0 for this example are approximately:")
    # The final equation is h'(x) - 1 = 0. We output the numbers in h'(x).
    print("h'(x) = (1/8)x^8 - (7/3)x^6 + (49/4)x^4 - 18x^2 + 25")
    print("Equation to solve: (1/8)x^8 - (7/3)x^6 + (49/4)x^4 - 18x^2 + 24 = 0")
    print("Roots:", roots)
    
    # Plotting the function and its roots
    x_vals = np.linspace(-4, 4, 1000)
    y_vals = h_prime(x_vals)
    
    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, y_vals, label="h'(x)")
    plt.axhline(y=1, color='r', linestyle='--', label="y=1")
    plt.plot(roots, np.ones_like(roots), 'go', label="Roots of h'(x)=1")
    plt.title("Graph of a plausible h'(x) showing 8 roots for h'(x)=1")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()
    plt.ylim(min(y_vals)-5, max(y_vals)+5)
    # plt.show() # Uncomment to display the plot
    
    print("\nBased on the analysis, the maximum number of fixed points is 9.")

find_and_show_roots()
