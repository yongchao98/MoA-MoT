import numpy as np
import matplotlib.pyplot as plt

def demonstrate_counterexample():
    """
    This function demonstrates a counterexample to the statement:
    'Any strictly convex function has a unique global minimizer'.

    The function f(x) = e^x is used as the counterexample.
    1. A function is strictly convex if its second derivative is always positive.
       For f(x) = e^x, the second derivative is f''(x) = e^x, which is > 0 for all real x.
       Therefore, f(x) = e^x is strictly convex.

    2. This code shows that f(x) = e^x does not have a global minimum on the set of real numbers.
    """
    print("Analyzing the strictly convex function f(x) = e^x:")
    print("--------------------------------------------------")
    print("The second derivative is f''(x) = e^x, which is always positive, proving strict convexity.")

    # Evaluate the function for increasingly negative values of x
    print("\nEvaluating the function as x approaches negative infinity:")
    x_values = [-1, -10, -50, -100]
    for x in x_values:
      # The final equation is f(x) = e^x. We output numbers for x and f(x).
      print(f"For x = {x}, the function value is f({x}) = {np.exp(x):.4e}")

    print("\nConclusion:")
    print("As x becomes more negative, f(x) gets closer and closer to 0 but never reaches it.")
    print("There is no single 'x' value that produces a minimum value for the function.")
    print("The infimum (greatest lower bound) is 0, but it is never attained.")
    print("This proves that the statement 'Any strictly convex function has a unique global minimizer' is FALSE.")

    # Plotting for visual confirmation
    x_plot = np.linspace(-5, 5, 400)
    y_plot = np.exp(x_plot)

    plt.figure(figsize=(8, 5))
    plt.plot(x_plot, y_plot, label='f(x) = e^x')
    plt.title('Plot of a Strictly Convex Function with No Global Minimum')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.grid(True)
    plt.legend()
    plt.axhline(0, color='red', linestyle='--', label='Infimum (y=0)')
    plt.ylim(bottom=-0.5, top=20)
    print("\nA plot will be generated to visualize this. Please close the plot window to continue.")
    plt.show()


demonstrate_counterexample()