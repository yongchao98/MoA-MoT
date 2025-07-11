import itertools

# Define the truth values using numeric representation for easier computation
T = 1.0   # True
G = 0.5   # Glut (both True and False)
F = 0.0   # False

# The set of all possible truth values
VALUES = [T, G, F]

# The set of designated values (values that count as "true" for validity)
DESIGNATED_VALUES = {T, G}

# --- Define the logical connectives for the 3-valued logic LP ---
def op_not(v):
    """Logical NOT (Negation)"""
    return 1.0 - v

def op_and(v1, v2):
    """Logical AND (Conjunction)"""
    return min(v1, v2)

def op_or(v1, v2):
    """Logical OR (Disjunction)"""
    return max(v1, v2)

def op_implies(v1, v2):
    """Logical IMPLIES (Material Implication), defined as not(v1) or v2"""
    return op_or(op_not(v1), v2)

def main():
    """
    Main function to evaluate the conclusion of argument L and check for tautology.
    Argument L: A |- (A & B) -> (B & A)
    """
    print("Evaluating the validity of argument L: A |- (A & B) -> (B & A)")
    print("The argument is valid if its conclusion is a tautology.")
    print("Let's check if the conclusion '(A & B) -> (B & A)' is a tautology in 3-valued logic.")
    print("\nTruth table for the conclusion:")
    print("-" * 65)
    print(f"{'A':<5} | {'B':<5} | {'A & B':<7} | {'B & A':<7} | {'(A & B) -> (B & A)':<20} | {'Is Designated?':<15}")
    print("-" * 65)

    is_tautology = True
    
    # Iterate through all possible truth value assignments for variables A and B
    for v_a in VALUES:
        for v_b in VALUES:
            # Evaluate the left side of the implication
            left_side = op_and(v_a, v_b)
            # Evaluate the right side of the implication
            right_side = op_and(v_b, v_a)
            # Evaluate the full conclusion formula
            result = op_implies(left_side, right_side)
            
            # Check if the result is a designated value
            is_designated = result in DESIGNATED_VALUES
            if not is_designated:
                is_tautology = False

            # Print the row of the truth table
            print(f"{v_a:<5.1f} | {v_b:<5.1f} | {left_side:<7.1f} | {right_side:<7.1f} | {result:<20.1f} | {str(is_designated):<15}")

    print("-" * 65)
    
    if is_tautology:
        print("\nThe conclusion is a tautology because its value is always designated (1.0 or 0.5).")
        print("An argument with a tautological conclusion is always valid.")
        print("Therefore, argument L is valid in system KG.")
    else:
        print("\nThe conclusion is not a tautology, so the argument's validity cannot be proven this way.")

    print("\nFinal formula checked (conclusion of L):")
    # Output the formula as requested by the prompt
    print("(A ∧ B) → (B ∧ A)")

if __name__ == '__main__':
    main()
