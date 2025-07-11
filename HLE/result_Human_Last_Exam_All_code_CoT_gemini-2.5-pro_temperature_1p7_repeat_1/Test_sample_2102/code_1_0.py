import sympy
import numpy as np

def solve():
    """
    Finds the smallest n where f(n) > 10 and computes n * ||W_n||_inf.
    """
    # Define the symbolic variable and the function g(x)
    x = sympy.Symbol('x')
    g = (2 / sympy.pi) * sympy.elliptic_k(x) * sympy.exp(x)

    # Set a tolerance for comparing floating-point roots to check for multiplicity
    TOLERANCE = 1e-9

    n = 1
    while True:
        # Step 1: Get the Taylor polynomial P_n(x) for g(x)
        # We need the series up to x^n, which has n+1 terms
        series = g.series(x, 0, n + 1).removeO()
        
        # Convert the sympy series to a sympy polynomial to extract coefficients easily
        poly = sympy.Poly(series, x)
        
        # Get coefficients from highest power to lowest for numpy.roots
        coeffs_sympy = poly.all_coeffs()
        coeffs_float = [float(c) for c in coeffs_sympy]

        # Step 2: Find the roots of the polynomial (eigenvalues of S_n)
        # np.roots finds the roots of a polynomial with coefficients from high to low power.
        if len(coeffs_float) <= 1: # Handles constant polynomial case
             roots = np.array([])
        else:
            roots = np.roots(coeffs_float)

        # Step 3: Calculate f(n)
        if roots.size == 0:
            f_n = 0
        else:
            f_n = np.sum(np.abs(roots)**3)

        # Step 4: Check if f(n) > 10
        if f_n > 10:
            print(f"Found the smallest n = {n}, where f(n) = {f_n:.4f} > 10.")
            
            # Step 5: Calculate ||W_n||_inf
            # This requires checking for repeated roots (multiplicities).
            
            # Group close roots together
            sorted_roots = sorted(roots, key=lambda z: (z.real, z.imag))
            unique_roots = []
            multiplicities = []
            
            if len(sorted_roots) > 0:
                current_root_cluster = [sorted_roots[0]]
                for i in range(1, len(sorted_roots)):
                    # If the next root is close to the first root of the current cluster
                    if np.abs(sorted_roots[i] - current_root_cluster[0]) < TOLERANCE:
                        current_root_cluster.append(sorted_roots[i])
                    else:
                        # Finalize the previous cluster
                        # The "unique root" is the average of the cluster
                        unique_roots.append(np.mean(current_root_cluster))
                        multiplicities.append(len(current_root_cluster))
                        
                        # Start a new cluster
                        current_root_cluster = [sorted_roots[i]]
                
                # Add the last cluster
                unique_roots.append(np.mean(current_root_cluster))
                multiplicities.append(len(current_root_cluster))

            # Calculate the infinity norm contribution from each unique eigenvalue block
            norm_contributions = []
            for i in range(len(unique_roots)):
                root = unique_roots[i]
                mult = multiplicities[i]
                
                # If multiplicity > 1, the norm of the block is |lambda| + 1.
                # Otherwise, it's just |lambda|.
                if mult > 1:
                    norm_contributions.append(np.abs(root) + 1)
                else:
                    norm_contributions.append(np.abs(root))
            
            # The infinity norm of the full Weyr matrix W_n is the maximum of the block norms.
            if not norm_contributions:
                inf_norm = 0
            else:
                inf_norm = max(norm_contributions)
                
            # Step 6: Final Calculation and Output
            result = n * inf_norm
            
            print(f"\nFor n = {n}, we calculate the infinity norm of W_n, ||W_n||_inf.")
            print(f"The unique eigenvalues have been identified along with their multiplicities.")
            print(f"The calculated value for ||W_n||_inf is: {inf_norm:.4f}")
            print("\nThe final requested value is n * ||W_n||_inf.")
            print(f"Equation: {n} * {inf_norm:.4f} = {result:.4f}")
            
            # Final answer in the required format
            print(f"\n<<<{result}>>>")
            break
        
        n += 1
        # Safety break to prevent infinite loops in case of unexpected behavior
        if n > 100:
            print("Loop limit reached, solution not found.")
            break

if __name__ == '__main__':
    solve()