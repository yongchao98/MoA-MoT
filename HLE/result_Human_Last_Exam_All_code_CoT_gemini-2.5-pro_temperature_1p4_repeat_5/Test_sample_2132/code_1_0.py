import math

def calculate_bottcher_complexity():
    """
    Calculates the Böttcher Molecular Complexity for the product of the
    Favorskii rearrangement of 2-chlorocyclohexanone.
    """
    # Step 1: Identify the product and its skeletal atoms.
    # The product is cyclopentanecarboxylic acid (C6H10O2).
    # The skeletal (non-hydrogen) atoms are 6 carbons and 2 oxygens.
    n_C = 6  # Number of Carbon atoms
    n_O = 2  # Number of Oxygen atoms

    # Step 2: Calculate the total number of skeletal atoms (N_skeletal).
    N_skeletal = n_C + n_O

    # Step 3: Define the Böttcher Molecular Complexity formula.
    # CI = N_skeletal*log2(N_skeletal) - SUM(n_i * log2(n_i))

    # Step 4: Calculate each term of the formula.
    log2_N = math.log2(N_skeletal)
    log2_C = math.log2(n_C)
    log2_O = math.log2(n_O)

    term1 = N_skeletal * log2_N
    term_C = n_C * log2_C
    term_O = n_O * log2_O

    # Step 5: Calculate the final complexity index.
    complexity = term1 - (term_C + term_O)

    # Step 6: Print the results and the equation.
    print("The product of the Favorskii rearrangement of 2-chlorocyclohexanone is cyclopentanecarboxylic acid (C6H10O2).")
    print("\nBöttcher Molecular Complexity Calculation:")
    print("Formula: CI = N_skeletal * log2(N_skeletal) - (n_C * log2(n_C) + n_O * log2(n_O))")
    print(f"\nValues:")
    print(f"  n_C (Carbon atoms) = {n_C}")
    print(f"  n_O (Oxygen atoms) = {n_O}")
    print(f"  N_skeletal (Total skeletal atoms) = {N_skeletal}")
    
    # Print the equation with all numbers substituted
    print("\nFinal Equation:")
    print(f"{N_skeletal} * log2({N_skeletal}) - ({n_C} * log2({n_C}) + {n_O} * log2({n_O}))")
    print(f"= {N_skeletal:.3f} * {log2_N:.3f} - ({n_C:.3f} * {log2_C:.3f} + {n_O:.3f} * {log2_O:.3f})")
    print(f"= {term1:.3f} - ({term_C:.3f} + {term_O:.3f})")
    print(f"= {term1:.3f} - {term_C + term_O:.3f}")
    print(f"= {complexity:.3f}")
    print(f"\nThe Böttcher Molecular Complexity is {complexity:.3f}")

if __name__ == "__main__":
    calculate_bottcher_complexity()
