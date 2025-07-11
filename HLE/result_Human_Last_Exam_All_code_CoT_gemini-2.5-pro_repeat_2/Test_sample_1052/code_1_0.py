def solve_boolean_puzzle():
    """
    This function presents the solution to the boolean formula puzzle.
    It shows the original Zhigalkin polynomial and the final derived Boolean formula.
    """

    # The given Zhigalkin polynomial in a readable format.
    polynomial_str = "c ⊕ d ⊕ (b∧c) ⊕ (a∧c) ⊕ (a∧d) ⊕ (a∧b∧c) ⊕ (a∧b∧d)"

    # The derived Boolean formula using the allowed set of operators.
    # The derivation steps are explained above.
    # Formula: (b AND d) OR ( (a NOR b) AND NOT(c XNOR d) )
    # Converted Formula: NOT(b -> NOT d) OR NOT( (a NOR b) -> (c XNOR d) )
    final_formula = "¬(b → ¬d) ∨ ¬((a ↓ b) → (c ↔ d))"

    print("The Zhigalkin polynomial is:")
    print(f"P = {polynomial_str}")
    print("\nThe corresponding Boolean formula using operators from the set (¬, ↑, ↓, ↔︎, →, ∨) is:")
    print(f"F(a,b,c,d) = {final_formula}")

solve_boolean_puzzle()