import numpy as np
import matplotlib.pyplot as plt

# Statement E is false.
# A strictly convex function is guaranteed to have AT MOST ONE global minimum, but it is not guaranteed to have one.
# For example, f(x) = e^x is strictly convex, but it never reaches a minimum value. It approaches an infimum of 0 but never attains it.
# This script demonstrates this counterexample.

# Generate x values
x = np.linspace(-5, 5, 400)
# Calculate y values for the strictly convex function f(x) = e^x
y = np.exp(x)

# Create the plot
fig, ax = plt.subplots()
ax.plot(x, y, label='f(x) = $e^x$')

# Add labels and title
ax.set_xlabel('x')
ax.set_ylabel('f(x)')
ax.set_title('Counterexample: A Strictly Convex Function with No Global Minimum')
ax.axhline(0, color='red', linestyle='--', label='Infimum (y=0)')
ax.legend()
ax.grid(True)
ax.set_ylim(bottom=-0.5)

# Add explanatory text
text_str = (
    "The function f(x) = e^x is strictly convex.\n"
    "However, it has no global minimum point.\n"
    "It approaches its greatest lower bound (infimum) of 0 as x approaches -∞,\n"
    "but it never actually reaches a value of 0 or less.\n"
    "Therefore, the statement 'Any strictly convex function has a unique global minimizer' is FALSE."
)
print(text_str)

# To avoid showing the plot in a real execution environment where it might not be supported,
# we save it to a file and print a confirmation message.
# In an interactive environment (like a Jupyter notebook), you could use plt.show().
try:
    plt.savefig("convex_function_plot.png")
    print("\nPlot saved to convex_function_plot.png")
except Exception as e:
    print(f"\nCould not save plot due to an error: {e}")
