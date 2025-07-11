import numpy as np

def T(x):
    """
    Implements the map T(x) = 1/sqrt(x) mod 1.
    """
    if x == 0:
        return 0
    # Add a check for robustness against invalid inputs
    if x < 0:
        raise ValueError("Input must be non-negative")
    inv_sqrt_x = 1.0 / np.sqrt(x)
    return inv_sqrt_x - np.floor(inv_sqrt_x)

def find_invariant_density():
    """
    Numerically finds the invariant density for the map T(x).
    """
    print("This script numerically estimates the invariant density rho(x) for the map T(x) = 1/sqrt(x) mod 1.")
    print("The density is assumed to have the functional form: rho(x) = C1/sqrt(x) + C2.")
    print("Generating a long orbit of the map to sample the measure...")

    # Step 1: Generate a long orbit
    N = 2000000  # Number of iterations for the orbit
    x = np.random.rand() # Start with a random point in (0, 1)

    # Discard the first few thousand points to let the orbit settle onto the attractor
    for _ in range(1000):
        x = T(x)
    
    # Use a list comprehension for efficient generation of the orbit
    orbit = np.array([x := T(x) for _ in range(N)])

    print("Creating a histogram and fitting the data...")
    # Step 2: Create a histogram to approximate the density
    num_bins = 500
    counts, bin_edges = np.histogram(orbit, bins=num_bins, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    
    # Exclude bins very close to x=0 to avoid singularity issues with the 1/sqrt(x) term during fitting
    mask = bin_centers > 1e-4
    bin_centers = bin_centers[mask]
    counts = counts[mask]
    
    # Step 3: Perform a linear least squares fit
    # We are fitting the model: counts = C1_un * (1/sqrt(bin_centers)) + C2_un * (1)
    y = counts
    x1 = 1.0 / np.sqrt(bin_centers)
    x2 = np.ones_like(bin_centers)
    M = np.vstack((x1, x2)).T # This is the design matrix
    
    try:
        # Solve the linear system for the unnormalized coefficients
        params, _, _, _ = np.linalg.lstsq(M, y, rcond=None)
        C1_un, C2_un = params
    except np.linalg.LinAlgError:
        print("Could not solve the least squares problem. Exiting.")
        return

    # Step 4: Normalize the found density function
    # The normalization constant is the integral of the unnormalized density from 0 to 1.
    # Integral[C1_un/sqrt(x) + C2_un] dx from 0 to 1 = [2*C1_un*sqrt(x) + C2_un*x] from 0 to 1
    # which evaluates to 2*C1_un + C2_un
    norm_const = 2 * C1_un + C2_un
    
    if norm_const == 0:
        print("Normalization constant is zero, cannot normalize. Exiting.")
        return
        
    C1_norm = C1_un / norm_const
    C2_norm = C2_un / norm_const

    # Step 5: Output the final equation
    print("\nThe estimated coefficients for the normalised density rho(x) = C1/sqrt(x) + C2 are:")
    print(f"C1 = {C1_norm:.4f}")
    print(f"C2 = {C2_norm:.4f}")
    print("\nThe final equation for the normalised density is:")
    print(f"rho(x) = {C1_norm:.4f} / (x**0.5) + {C2_norm:.4f}")

# Run the function to find the density
find_invariant_density()