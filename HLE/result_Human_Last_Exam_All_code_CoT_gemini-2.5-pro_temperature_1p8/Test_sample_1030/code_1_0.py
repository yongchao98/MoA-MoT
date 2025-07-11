def check_kg_logic():
    """
    This script defines and evaluates logic expressions in the 3-valued logic KG.
    KG has truth values {True, Glut, False} and its designated values are {True, Glut}.
    The script checks the validity of argument L and prints the evaluation trace.
    """
    
    # Define truth values: T=2 (True), G=1 (Glut), F=0 (False)
    vals = [2, 1, 0]

    def v_neg(v):
        """Computes negation."""
        if v == 2: return 0  # ¬T = F
        if v == 0: return 2  # ¬F = T
        return 1           # ¬G = G

    def v_and(v1, v2):
        """Computes conjunction (min value)."""
        return min(v1, v2)

    def v_or(v1, v2):
        """Computes disjunction (max value)."""
        return max(v1, v2)

    def v_imp(v1, v2):
        """Computes implication A → B as ¬A ∨ B."""
        return v_or(v_neg(v1), v2)

    def is_designated(v):
        """Checks if a value is designated (T or G)."""
        return v >= 1

    def to_char(v):
        """Converts numeric value to character representation."""
        if v == 2: return 'T'
        if v == 1: return 'G'
        return 'F'

    # Check argument L: A ⊢ (A ∧ B) → (B ∧ A)
    # An argument is valid if there's no case where all premises are designated
    # and the conclusion is not.

    print("Checking validity of argument L: A ⊢ (A ∧ B) → (B ∧ A)")
    print("Designated truth values are T(2) and G(1).")
    print("-" * 75)
    print("Inputs | Premise |       Conclusion Calculation      | Result of Check")
    print(" A | B |   A   | A ∧ B | B ∧ A | (A ∧ B) → (B ∧ A) | (Premise Desig.) → (Concl. Desig.)?")
    print("-" * 75)

    counter_example_found = False
    for a_val in vals:
        for b_val in vals:
            # Evaluate the premise
            premise_val = a_val
            
            # Evaluate the conclusion
            term1_val = v_and(a_val, b_val)
            term2_val = v_and(b_val, a_val)
            conclusion_val = v_imp(term1_val, term2_val)

            # Check for invalidity
            validity_char = '✓'
            if is_designated(premise_val) and not is_designated(conclusion_val):
                counter_example_found = True
                validity_char = 'X <-- COUNTER-EXAMPLE'

            # Print each number (as a character) in the final equation
            print(f" {to_char(a_val)} | {to_char(b_val)} |   {to_char(premise_val)}   |   {to_char(term1_val)}   |   {to_char(term2_val)}   |         {to_char(conclusion_val)}         | {validity_char}")

    print("-" * 75)
    if not counter_example_found:
        print("Result: No counter-examples found. The argument is valid.")
        print("\nNote: The conclusion column shows only designated values (T or G).")
        print("This means the conclusion is a tautology in KG. An argument with a tautological conclusion is always valid.")
    else:
        print("Result: A counter-example was found. The argument is invalid.")

# Execute the check
check_kg_logic()