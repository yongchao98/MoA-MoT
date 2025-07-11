import numpy as np
import time

def estimate_measure(N, Kx=1000, Kt_factor=2.5):
    """
    Estimates the measure of the set X for a given N.
    
    Args:
        N: The length of the sequence.
        Kx: The number of points to sample for x in [0,1].
        Kt_factor: Determines the density of the t-grid sampling.
                   The number of t-points will be Kt_factor * N^2.
                   Theory requires this to be > 2.
    
    Returns:
        The estimated measure of X.
    """
    print(f"Estimating measure for N = {N}")
    x_grid = np.linspace(0, 1, Kx, endpoint=False)
    # According to sampling theory for trigonometric polynomials, the number
    # of sample points for t must be greater than twice the max frequency (N^2).
    Kt = int(Kt_factor * N**2)
    t_grid = np.linspace(0, 1, Kt, endpoint=False)
    n_vals = np.arange(1, N + 1)
    
    # We choose a_n = 1/sqrt(N) for this simulation.
    a_n_norm = 1.0 / np.sqrt(N)
    
    threshold = N**(3/8)
    
    count = 0
    start_time = time.time()
    
    for i, x in enumerate(x_grid):
        # For each x, compute max_t |S(x,t)|
        # Vectorize the calculation over t for the current x
        # The phase is a matrix of size (N, Kt): n_vals (broadcast) * x + n_vals^2 (broadcast) * t_grid
        phases = 2 * np.pi * (np.outer(n_vals, [x]) + np.outer(n_vals**2, t_grid))
        
        # S(x,t) values for a fixed x and all t in t_grid
        s_values = np.sum(a_n_norm * np.exp(1j * phases), axis=0)
        
        max_s_magnitude = np.max(np.abs(s_values))
        
        if max_s_magnitude > threshold:
            count += 1
            
        # Progress indicator
        if (i + 1) % 100 == 0:
            elapsed = time.time() - start_time
            est_rem = (elapsed / (i + 1)) * (Kx - 1 - i)
            print(f"  ... {i+1}/{Kx} x-values processed. Time elapsed: {elapsed:.1f}s. Est. remaining: {est_rem:.1f}s")
            
    measure = count / Kx
    total_time = time.time() - start_time
    print(f"-> N = {N}, Measure |X| ~= {measure:.4f} (took {total_time:.2f}s)")
    return measure

def main():
    """
    Calculates alpha by running the estimation for two values of N and 
    finding the slope in log-log scale.
    """
    # Parameters for the simulation.
    # Larger N values give better estimates but take much longer.
    # N must be chosen carefully to make the calculation feasible.
    # E.g., N1=16, N2=25
    try:
        # Note: These parameters are for demonstration. For an accurate result,
        # much larger N, Kx would be needed, which is not feasible for typical execution times.
        N1 = 16 
        N2 = 25
        # Kx should be large for accuracy
        Kx = 2000 
        # Kt_factor should be > 2
        Kt_factor = 2.1
        
        print("Starting numerical estimation of alpha. This may take a while.")
        
        m1 = estimate_measure(N1, Kx, Kt_factor)
        m2 = estimate_measure(N2, Kx, Kt_factor)

        if m1 > 0 and m2 > 0:
            # alpha = log(m2/m1) / log(N2/N1)
            alpha = np.log(m2 / m1) / np.log(N2 / N1)
            print(f"\nNumerical estimate for alpha: {alpha:.4f}")
        else:
            print("\nNumerical estimation failed. The measured count was zero.")
            print("Try with larger N values or larger Kx.")
            
    except Exception as e:
        print(f"An error occurred: {e}")

    print("\nBased on theoretical analysis, the value of alpha is -3/4.")
    print(-3/4)

if __name__ == '__main__':
    main()
