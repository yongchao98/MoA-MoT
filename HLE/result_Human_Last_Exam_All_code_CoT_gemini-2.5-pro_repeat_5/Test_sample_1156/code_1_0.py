import numpy as np

def solve():
    """
    Calculates and prints a numerical approximation of the normalized 
    invariant density for the map T(x) = 1/sqrt(x) mod 1.
    """

    # The map T(x) = 1/sqrt(x) mod 1
    def T(x):
      if x == 0:
        return 0.0
      return (1.0 / np.sqrt(x)) % 1.0

    # Simulation parameters
    N_ITERATIONS = 10**7  # Number of iterations for the simulation
    N_BINS = 50           # Number of bins for the histogram
    
    # We start with a random point to generate a typical trajectory.
    # We run a few iterations first to let the trajectory settle onto the attractor.
    x = np.random.rand()
    for _ in range(100):
        x = T(x)

    # Generate the trajectory and store the points
    trajectory = np.zeros(N_ITERATIONS)
    for i in range(N_ITERATIONS):
      x = T(x)
      trajectory[i] = x

    # Compute the histogram of the trajectory points
    hist, bin_edges = np.histogram(trajectory, bins=N_BINS, range=(0.0, 1.0))
    
    # Normalize the histogram to get the density
    # The density rho(x_i) is approximated by the count in the bin
    # divided by the total number of points and the bin width.
    bin_width = bin_edges[1] - bin_edges[0]
    density = hist / (N_ITERATIONS * bin_width)
    
    # Get the center of each bin for printing
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    # The "final equation" is interpreted as the set of values that
    # define our approximate density function, rho(x_i) = y_i.
    # We print each of these numbers.
    print("Approximated density rho(x) at the center of each bin:")
    print("-----------------------------------------------------")
    print("    x     |   rho(x)   ")
    print("-----------------------------------------------------")
    for center, d_val in zip(bin_centers, density):
        # We output the numbers defining the approximate density function
        # in the format "rho(center) = density_value".
        print(f"rho({center:.4f}) = {d_val:.6f}")

solve()