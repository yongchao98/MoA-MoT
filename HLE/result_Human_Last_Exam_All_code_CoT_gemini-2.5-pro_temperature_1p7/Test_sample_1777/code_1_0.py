import sys

# Define a function to solve the problem
def calculate_rydberg_energy():
    """
    Calculates the Rydberg energy for the n=3 state of an exciton
    in a 2D semiconductor.
    """
    # --- Given Parameters ---
    # Band gap in eV
    E_g = 3.0
    # 1s exciton resonance peak in eV
    E_1s = 1.0
    # Principal quantum number for the target state
    n = 3

    # --- Step 1: Calculate the ground state (1s) binding energy ---
    # The exciton energy is E_n = E_g - E_bn, so the binding energy is E_bn = E_g - E_n.
    E_b1 = E_g - E_1s

    print("Step 1: Calculate the ground state (1s) binding energy (E_b1).")
    print(f"The formula is: E_b1 = E_g - E_1s")
    print(f"E_b1 = {E_g} eV - {E_1s} eV = {E_b1} eV")
    print("-" * 30)

    # --- Step 2: Calculate the Rydberg energy for n=3 ---
    # For a 2D system, the binding energy formula is E_bn = E_b1 / (4 * (n - 0.5)^2)
    numerator = E_b1
    n_minus_half = n - 0.5
    denominator_part = (n - 0.5)**2
    denominator_full = 4 * denominator_part
    E_bn = numerator / denominator_full

    print(f"Step 2: Calculate the Rydberg energy (binding energy) for n = {n} (E_b{n}).")
    print("The formula for a 2D exciton is: E_bn = E_b1 / (4 * (n - 0.5)^2)")
    print("\nThe calculation for the final equation is as follows:")
    print(f"E_b{n} = {numerator} / (4 * ({n} - 0.5)^2)")
    print(f"E_b{n} = {numerator} / (4 * {n_minus_half}^2)")
    print(f"E_b{n} = {numerator} / (4 * {denominator_part})")
    print(f"E_b{n} = {numerator} / {denominator_full}")
    print(f"E_b{n} = {E_bn:.3f} eV")

    # This is to be captured by the system for the final answer format
    sys.stdout.write(f"\n<<<{E_bn:.2f}>>>")

# Execute the function
calculate_rydberg_energy()