import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    This program calculates and visualizes an approximation of the 
    invariant density for the map T(x) = 1/sqrt(x) mod 1.
    """
    
    # Define the map T(x)
    def T(x):
      if x == 0:
        return 0
      # np.modf separates the fractional and integer parts. We need the fractional part.
      return np.modf(x**(-0.5))[0]

    # --- Simulation ---
    # Set a random seed for reproducibility
    np.random.seed(0)
    # Start with a random initial value
    x = np.random.rand()
    # Number of iterations for the simulation
    num_iterations = 2000000

    # Create an array to store the trajectory
    trajectory = np.zeros(num_iterations)
    # Run the simulation
    for i in range(num_iterations):
      x = T(x)
      trajectory[i] = x

    # --- Plotting ---
    plt.figure(figsize=(12, 7))
    
    # Plot the histogram of the trajectory
    plt.hist(trajectory, bins=250, density=True, color='lightblue', label='Numerical Simulation (Histogram)')

    # --- Analytical Approximation ---
    # Define the proposed approximate density function rho(x) = (4/11)*x^(-1/2) + 3/11
    coeff_a = 4.0 / 11.0
    coeff_b = 3.0 / 11.0
    
    def rho(x):
      # We avoid x=0 directly to prevent division by zero warning
      # The integral is finite, so this is just for plotting
      safe_x = np.maximum(x, 1e-9)
      return coeff_a * safe_x**(-0.5) + coeff_b

    # Generate x values for the plot of the analytical approximation
    x_vals = np.linspace(0.001, 1, 500)
    y_vals = rho(x_vals)

    # Plot the analytical approximation
    plt.plot(x_vals, y_vals, 'r-', linewidth=2, label=r'Approximation: $\rho(x) = \frac{4}{11}x^{-1/2} + \frac{3}{11}$')

    plt.title(r'Invariant Density for $T(x) = x^{-1/2} \mod 1$')
    plt.xlabel('x')
    plt.ylabel(r'Density $\rho(x)$')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.ylim(0, 3.5) # Limit y-axis for better visibility
    plt.show()

    # --- Outputting the numbers in the final equation ---
    print("The approximate normalised density is of the form: A * x^(-1/2) + B")
    print("The coefficients are:")
    print(f"A = {coeff_a}")
    print(f"B = {coeff_b}")
    print("\nFinal Equation: rho(x) = (4/11) * x^(-1/2) + (3/11)")
    print("Each number in the final equation is:")
    print(4)
    print(11)
    print(-1)
    print(2)
    print(3)
    print(11)

if __name__ == '__main__':
    main()
