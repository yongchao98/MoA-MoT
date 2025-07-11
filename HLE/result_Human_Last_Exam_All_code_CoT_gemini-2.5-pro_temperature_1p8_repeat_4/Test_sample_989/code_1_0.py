import numpy as np
import matplotlib.pyplot as plt

def demonstrate_false_statement_E():
    """
    This function demonstrates that statement E is false by plotting a counterexample.
    Statement E: Any strictly convex function has a unique global minimizer.
    This is false because the existence of a minimizer is not guaranteed.
    """

    # Define a strictly convex function: f(x) = e^x
    x = np.linspace(-10, 3, 400)
    y = np.exp(x)

    # Explanation
    print("Analyzing Statement E: 'Any strictly convex function has a unique global minimizer.'")
    print("\nThis statement is FALSE.")
    print("A function can be strictly convex but not have a global minimum.")
    print("A classic counterexample is the exponential function, f(x) = e^x.")
    print("\nThe plot below shows f(x) = e^x. As you can see, the function is convex (curves upwards).")
    print("As x approaches negative infinity, the function value approaches 0 but never reaches it.")
    print("Therefore, this function does not have a global minimum, which disproves the statement.")

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, label='f(x) = e^x')

    # Add plot elements for clarity
    plt.title('Counterexample for Statement E')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.axhline(0, color='red', linestyle='--', label='Infimum (y=0)')
    plt.grid(True)
    plt.legend()
    plt.ylim(-1, 15)
    plt.show()

if __name__ == '__main__':
    demonstrate_false_statement_E()