import numpy as np
import matplotlib.pyplot as plt

def demonstrate_convex_no_minimum():
    """
    This function demonstrates that a strictly convex function does not necessarily have a global minimum.
    It plots f(x) = e^x, which is a strictly convex function.
    """
    print("Analyzing the statement: 'Any strictly convex function has a unique global minimizer.'")
    print("\nThis statement is mathematically false.")
    print("A counterexample is the function f(x) = e^x. This function is strictly convex.")
    print("However, it does not have a global minimum on the set of real numbers.")
    print("As x approaches negative infinity, f(x) approaches 0, but it never reaches 0.")
    print("The plot below visualizes this function.")

    # Generate x values
    x = np.linspace(-10, 5, 400)
    # Calculate y values for f(x) = e^x
    y = np.exp(x)

    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, label='f(x) = e^x')
    plt.title('Plot of a Strictly Convex Function with No Minimum')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.grid(True)
    plt.legend()
    plt.axhline(0, color='red', linestyle='--', label='Infimum y=0 (never reached)')
    plt.legend()
    
    # Add text to explain the plot
    plt.text(-9.5, 20000, "The function f(x)=e^x is strictly convex.", fontsize=10)
    plt.text(-9.5, 12000, "Its value approaches 0 as x -> -inf, but never reaches a minimum.", fontsize=10)

    print("\nDisplaying plot...")
    plt.show()

if __name__ == '__main__':
    demonstrate_convex_no_minimum()