def generate_formulas():
    """
    This function prints the definitions for the function f(w) and the formula C
    that encode the equipartitioning problem in linear logic.
    """

    # Define the components of the formulas using only allowed elements.
    # We build a base "atomic" formula, A, that is not provable from an empty context.
    atom_A_str = "((1 ⊸ ⊥) ⊸ ⊥)"

    # We build a second "atomic" formula, S, distinct from A.
    atom_S_str = f"({atom_A_str} ⊸ ⊥)"

    # Explain the function f(w)
    print("--- 1. The Function f(w) ---")
    print("The function f(w) maps a natural number w to a linear logic formula.")
    print("Let a base 'atomic' formula A be defined as:")
    print(f"  A = {atom_A_str}")
    print("\nf(w) is defined as the w-th tensor power of A (for w > 0). f(0) is the multiplicative unit 1.")
    print("  f(w) = A ⊗ A ⊗ ... ⊗ A  (w times)")
    # We use ^ to denote tensor power for readability.
    print("  Symbolically: f(w) = A^w")
    print(f"  Example for w=4: f(4) = ({atom_A_str})^(4)")
    print("-" * 40)

    # Explain the formula C(W, m, b)
    print("\n--- 2. The Formula C(W, m, b) ---")
    print("The formula C defines the goal of the sequent.")
    print("Let a 'success token' formula S be defined as S = A ⊸ ⊥:")
    print(f"  S = {atom_S_str}")
    print("\nFirst, we define a 'bucket' formula B, representing the requirement for a single partition to sum to b:")
    print("  B = ((A^b) ⊸ S) ⊸ S")
    print(f"  Example for b=10: B = (({atom_A_str})^(10) ⊸ {atom_S_str}) ⊸ {atom_S_str}")
    print("\nThe full goal formula C is then m tensor products of this bucket formula B:")
    print("  C(W, m, b) = B ⊗ B ⊗ ... ⊗ B  (m times)")
    print("  Symbolically: C(W, m, b) = B^m")
    print(f"  Example for m=3, b=10: C = ( (({atom_A_str})^(10) ⊸ {atom_S_str}) ⊸ {atom_S_str} )^(3)")
    print("-" * 40)

    # State the final sequent
    print("\n--- 3. The Final Sequent ---")
    print("The equipartitioning problem EP(W, m, b) is true if and only if the following sequent is derivable in linear logic:")
    print("\n  { f(w) | w ∈ W } ⊢ C(W, m, b)\n")

# Execute the function to print the solution.
generate_formulas()