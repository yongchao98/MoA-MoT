import numpy as np

def H_exact(s, s_prime, R):
    """
    Calculates the exact chord length (generating function) for a circle.
    """
    x_s, y_s = R * np.cos(s / R), R * np.sin(s / R)
    x_sp, y_sp = R * np.cos(s_prime / R), R * np.sin(s_prime / R)
    return np.sqrt((x_sp - x_s)**2 + (y_sp - y_s)**2)

def H_asymptotic_and_components(s, s_prime, R):
    """
    Calculates the asymptotic approximation of the generating function H(s, s').
    It also returns the components of the formula for inspection.
    """
    ds = abs(s_prime - s)
    
    # For a circle, curvature is constant: kappa = 1/R
    kappa = 1.0 / R
    
    # Leading term
    term1 = ds
    
    # Leading-order correction term
    coeff_num = 1
    coeff_den = 24
    correction_term = (coeff_num / coeff_den) * (kappa**2) * (ds**3)
    
    # Asymptotic approximation
    H_approx = term1 - correction_term
    
    return H_approx, term1, coeff_num, coeff_den, kappa, ds

def main():
    """
    Main function to perform the analysis and print results.
    """
    R = 10.0  # Radius of the circular billiard table
    s = 0.0   # Starting point on the boundary

    print(f"Analyzing the generating function H(s, s') for a circle of radius R = {R}")
    print(f"The constant curvature is κ = 1/R = {1.0/R:.2f}")
    print("The theoretical asymptotic formula is: H(s, s') ≈ |s'-s| - (1/24)κ²|s'-s|³\n")
    print("-" * 75)

    # Test for progressively smaller arc-length separations (ds)
    for power in range(1, 5):
        ds = 10.0**(-power)
        s_prime = s + ds

        # Get exact and asymptotic values
        exact_val = H_exact(s, s_prime, R)
        asymptotic_val, t1, c_num, c_den, k_val, ds_val = H_asymptotic_and_components(s, s_prime, R)
        
        # Calculate error
        error = abs(exact_val - asymptotic_val)
        
        print(f"For arc-length separation |s'-s| = {ds_val:.1e}:")
        
        # This section outputs each number in the final equation as requested
        print(f"  Asymptotic Equation: H ≈ {t1:.9f} - ({c_num}/{c_den}) * ({k_val:.2f}**2) * ({ds_val:.1e}**3)")
        
        print(f"  -> H_asymptotic = {asymptotic_val:.9f}")
        print(f"  -> H_exact      = {exact_val:.9f}")
        print(f"  -> Error (Exact - Approx) = {error:.3e}")
        # The ratio of errors should be approx (1e-1)^5 / (1e-2)^5 = 1e5 etc.
        # This shows the error scales with ds^5
        if 'last_error' in locals() and last_error > 1e-15:
            print(f"  -> Error Reduction Factor ≈ {last_error/error:.1f} (Expected ≈ 100000.0 for O(ds^5))")
        print("-" * 75)

        last_error = error

if __name__ == "__main__":
    main()
