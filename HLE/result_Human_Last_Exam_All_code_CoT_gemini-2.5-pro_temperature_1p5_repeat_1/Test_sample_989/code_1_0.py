import numpy as np
import matplotlib.pyplot as plt

def demonstrate_false_statement_E():
    """
    This function demonstrates why statement E is false by plotting a counterexample.
    Statement E: "Any strictly convex function has a unique global minimizer."

    We will use the function f(x) = e^x.
    A function is strictly convex if its second derivative is strictly positive.
    The second derivative of e^x is e^x, which is always positive for any real x.
    Therefore, f(x) = e^x is a strictly convex function.

    This code will plot the function and show that it does not have a global minimizer.
    """

    # Generate x values
    x = np.linspace(-10, 5, 400)
    # Calculate y values for f(x) = e^x
    y = np.exp(x)

    print("--- Analysis of Statement E ---")
    print("Statement E: 'Any strictly convex function has a unique global minimizer.'")
    print("\nThis statement is FALSE.")
    print("Let's consider the function f(x) = e^x, which is strictly convex because its second derivative (e^x) is always positive.")
    print("The plot below shows that as x approaches negative infinity, f(x) approaches 0 but never reaches it.")
    print("Therefore, this function has no point 'x' that corresponds to a minimum value, so it has no global minimizer.")
    print("\nConclusion: Not all strictly convex functions have a global minimizer.")

    # Due to the environment, we can't show the plot directly with plt.show().
    # The printed explanation and the code logic serve as the main output.
    # To see the plot, you can run this code in a local Python environment
    # and uncomment the plt.show() line.

    # plt.figure(figsize=(10, 6))
    # plt.plot(x, y, label='f(x) = e^x (Strictly Convex)')
    # plt.axhline(0, color='red', linestyle='--', label='Infimum (y=0), but never reached')
    # plt.title('Counterexample for Statement E')
    # plt.xlabel('x')
    # plt.ylabel('f(x)')
    # plt.grid(True)
    # plt.legend()
    # plt.show()


demonstrate_false_statement_E()