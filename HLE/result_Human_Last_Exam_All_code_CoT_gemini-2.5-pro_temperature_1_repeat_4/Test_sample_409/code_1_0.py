import sys

def solve_cohomology_question():
    """
    This function explains the reasoning for determining the minimal cohomology degree
    for non-trivial extensions and obstructions in semi-abelian categories.
    """

    print("Analyzing the significance of cohomology degrees H^n(B, M):")
    print("=" * 60)

    # Degree 0
    degree_0 = 0
    print(f"For degree n = {degree_0}:")
    print(f"H^{degree_0}(B, M) represents the invariants of M under the action of B.")
    print("This is fundamental but does not classify extensions or obstructions.\n")

    # Degree 1
    degree_1 = 1
    print(f"For degree n = {degree_1}:")
    print(f"H^{degree_1}(B, M) classifies derivations and is related to split extensions.")
    print("While it is the first degree related to extensions, it does not capture the problem of general (non-split) extensions, nor is it the primary locus for obstructions.\n")

    # Degree 2
    degree_2 = 2
    print(f"For degree n = {degree_2}:")
    print(f"H^{degree_2}(B, M) is the first degree where both concepts become significant:")
    print(f"  1. Extensions: It classifies the set of (non-split) extensions of the object B by the B-module M.")
    print(f"  2. Obstructions: It is the home of the primary obstruction class in many lifting and deformation problems.")
    print("This makes it the minimal degree for the phenomena described.\n")

    # Degree 3
    degree_3 = 3
    print(f"For degree n = {degree_3}:")
    print(f"H^{degree_3}(B, M) and higher groups classify more complex structures and higher-order obstructions, but are not the minimal degree for the initial problem.\n")

    print("=" * 60)
    print("Conclusion: The minimal cohomology degree where non-trivial extensions and obstructions become significant is 2.")

    # Fulfilling the constraint to output each number in the final equation.
    final_equation_lhs = "Minimal Degree"
    final_equation_rhs = 2
    print(f"\nFinal Equation: {final_equation_lhs} = {final_equation_rhs}")

solve_cohomology_question()