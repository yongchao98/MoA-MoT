import numpy as np
import matplotlib.pyplot as plt

def explain_and_plot_counterexample():
    """
    Explains and plots a counterexample to the statement:
    'Any strictly convex function has a unique global minimizer.'
    """

    # --- Explanation ---
    print("Analyzing the statement: 'E. Any strictly convex function has a unique global minimizer.'")
    print("-" * 70)
    print("This statement is FALSE.")
    print("\nA function f(x) is strictly convex if its second derivative, f''(x), is greater than 0.")
    print("While it is true that IF a strictly convex function has a global minimum, it must be unique,")
    print("it is NOT guaranteed that a global minimum exists.\n")

    print("A classic counterexample is the function f(x) = e^x (the exponential function).")
    print("1. First derivative: f'(x) = e^x")
    print("2. Second derivative: f''(x) = e^x")
    print("\nSince f''(x) = e^x is always > 0 for any real number x, the function is strictly convex.")
    print("However, as x approaches negative infinity, f(x) approaches 0 but never reaches it.")
    print("The function does not have a single point that is a global minimum.")
    print("Therefore, we have found a strictly convex function that does not have a global minimizer,")
    print("which proves the original statement is false.\n")
    print("The plot below visualizes this function.")

    # --- Plotting ---
    x = np.linspace(-10, 3, 400)
    y = np.exp(x)

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, label='f(x) = e^x')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.title('Counterexample: A Strictly Convex Function with No Global Minimum')
    plt.xlabel('x')
    plt.ylabel('f(x) = e^x')
    plt.axhline(y=0, color='r', linestyle='--', label='Infimum (y=0, never reached)')
    plt.legend()
    plt.ylim(-0.5, 20)
    plt.text(-9.5, 10,
             "f''(x) = e^x > 0, so the function is strictly convex.",
             fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
    plt.text(-9.5, 7,
             "The function approaches 0 as x -> -inf,\nbut never attains a minimum value.",
             fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

    # In a real environment, plt.show() would be called.
    # For this exercise, we save it to a file.
    plot_filename = "convex_function_plot.png"
    plt.savefig(plot_filename)
    print(f"\nA plot illustrating this concept has been saved as '{plot_filename}'.")


if __name__ == '__main__':
    explain_and_plot_counterexample()
