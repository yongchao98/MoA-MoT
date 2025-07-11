import numpy as np

def find_invariant_density():
    """
    Numerically calculates the invariant density for the map T(x) = 1/sqrt(x) mod 1.
    The method uses iteration of the Perron-Frobenius operator.
    The final density is approximated by a polynomial and its equation is printed.
    """
    
    # --- 1. Set up the numerical parameters ---
    N = 2000      # Number of grid points for the density function
    K = 500       # Number of terms to include in the sum over preimages
    n_iter = 150  # Number of iterations
    tolerance = 1e-7 # Convergence tolerance
    poly_degree = 8  # Degree of the polynomial for fitting the result
    
    # --- 2. Initialize the density function ---
    grid = np.linspace(0, 1, N) # The grid on [0,1)
    rho = np.ones(N)            # Start with a uniform density rho_0(x) = 1
    
    # Normalize the initial density
    integral = np.trapz(rho, grid)
    rho = rho / integral
    
    print("Starting iterative calculation of the invariant density...")
    print("This may take a moment...")
    
    # --- 3. Iterate the Perron-Frobenius operator ---
    for it in range(n_iter):
        rho_new = np.zeros(N)
        for i, y in enumerate(grid):
            
            # The sum is over k=1, 2, ...
            k_vals = np.arange(1, K + 1)
            
            # Find the preimages x_k for the current y
            # We need to handle y=0 for the first term where k=0 for some maps, but here k>=1.
            # to avoid division by zero if y=0, we use a small epsilon, though not strictly needed here as k starts at 1
            y_safe = y if y > 0 else 1e-12
            preimages_x = 1.0 / (y_safe + k_vals)**2
            
            # Evaluate rho at these preimages using linear interpolation
            # np.interp expects the x-coordinates of the data points (our grid) to be increasing.
            # The preimages are decreasing in k, but we can pass them as is.
            rho_vals_at_preimages = np.interp(preimages_x, grid, rho)
            
            # The denominators in the Perron-Frobenius sum
            denominators = (y_safe + k_vals)**3
            
            # Calculate the sum
            s = 2 * np.sum(rho_vals_at_preimages / denominators)
            rho_new[i] = s

        # Normalize the new density
        integral = np.trapz(rho_new, grid)
        if integral > 1e-9:
          rho_new = rho_new / integral
    
        # Check for convergence
        diff = np.linalg.norm(rho_new - rho, ord=np.inf)
        if it % 10 == 0:
            print(f"Iteration {it+1}/{n_iter}, Change in density: {diff:.4g}")

        rho = rho_new
        if diff < tolerance:
            print(f"\nConverged after {it+1} iterations.")
            break
            
    if it == n_iter-1:
        print(f"\nReached maximum iterations ({n_iter}) without full convergence.")

    # --- 4. Fit the result to a polynomial ---
    print("\nFitting the numerical density to a polynomial...")
    coeffs = np.polyfit(grid, rho, poly_degree)
    
    # --- 5. Print the final equation ---
    equation_parts = []
    for i, c in enumerate(coeffs):
        power = poly_degree - i
        # Format the term based on its power
        if abs(c) < 1e-6: # neglect very small coefficients
            continue
            
        sign = "-" if c < 0 else "+"
        val = abs(c)

        if power > 1:
            term = f" {sign} {val:.4f}*x^{power}"
        elif power == 1:
            term = f" {sign} {val:.4f}*x"
        else: # power is 0
            term = f" {sign} {val:.4f}"
        equation_parts.append(term)
        
    final_equation = "rho(x) approx = " + "".join(equation_parts).lstrip(" +")
    print("\nAn approximate analytical expression for the normalised density is:")
    print(final_equation)
    
    # Print the peak value for the final answer
    peak_value = np.max(rho)
    return peak_value

if __name__ == '__main__':
    peak_density = find_invariant_density()
    # The final answer format requires a single value.
    # We choose the peak of the density function as a representative value.
    # To conform to the output format, we wrap the final numeric result in '<<<' and '>>>'.
    # This print will not appear if imported, only when run as a script.
    print(f"\nPeak value of the density function: {peak_density}")
    print(f"\n<<< {peak_density:.4f} >>>")
