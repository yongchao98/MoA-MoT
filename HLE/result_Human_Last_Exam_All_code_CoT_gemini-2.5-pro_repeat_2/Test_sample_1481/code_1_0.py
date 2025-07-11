import math

def demonstrate_asymptotic_H():
    """
    Demonstrates the asymptotic formula for the billiard generating function H(s, s')
    by comparing it to the exact value for a circular billiard table.
    """
    # --- Parameters ---
    # Radius of the circular billiard table
    R = 5.0
    # Arc-length separation (s' - s)
    ds = 0.5

    # --- Theoretical Values ---
    # For a circle, the curvature kappa is 1/R
    kappa = 1.0 / R
    
    # The asymptotic formula is H ≈ |s'-s| - (κ(s)²/24) * |s'-s|³
    # Let's calculate each component
    term1 = ds
    kappa_squared_over_24 = (kappa**2) / 24.0
    ds_cubed = ds**3
    
    # Calculate the final approximation
    H_approx = term1 - kappa_squared_over_24 * ds_cubed

    # --- Exact Calculation for a Circle ---
    # The exact chord length H is 2*R*sin(|s'-s| / (2*R))
    H_exact = 2 * R * math.sin(ds / (2 * R))

    # --- Output Results ---
    print("--- Asymptotic Analysis of Billiard Generating Function H(s, s') ---")
    print(f"\nAnalysis for a circle of Radius R = {R}")
    print(f"The curvature κ = 1/R = {kappa:.4f}")
    print(f"Separation |s' - s| = {ds}")
    
    print("\nAsymptotic Formula:")
    print("H(s, s') ≈ |s' - s| - (κ² / 24) * |s' - s|³")
    
    print("\nPlugging in the numbers:")
    # The user request: "Remember in the final code you still need to output each number in the final equation!"
    print(f"H ≈ {term1} - ({kappa**2:.4f} / 24) * {ds**3:.4f}")
    print(f"H ≈ {term1} - {kappa_squared_over_24:.6f} * {ds_cubed:.4f}")
    print(f"H ≈ {term1} - {kappa_squared_over_24 * ds_cubed:.6f}")
    
    print("\n--- Final Values ---")
    print(f"Asymptotic H ≈ {H_approx:.12f}")
    print(f"Exact H      = {H_exact:.12f}")
    print(f"Approximation Error = {abs(H_exact - H_approx):.4e}")

if __name__ == '__main__':
    demonstrate_asymptotic_H()