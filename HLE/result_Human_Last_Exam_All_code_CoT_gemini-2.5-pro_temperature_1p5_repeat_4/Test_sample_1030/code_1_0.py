#
# Plan:
# 1. Define the 3-valued logic system (KG) connectives.
#    - We will represent False as 0, Glut as 1, and True as 2.
#    - Designated values are {1, 2} (Glut, True).
# 2. Define the logical connectives NOT, AND, OR based on the rules of
#    3-valued logic.
# 3. The validity of argument L, A |- (A & B) -> (B & A), under strict
#    implication hinges on the inner formula P = `¬(A ∧ B) ∨ (B ∧ A)` being a
#    tautology (always having a designated value).
# 4. We will iterate through all 9 possible combinations of truth values for A and B.
# 5. For each combination, we calculate the truth value of P and verify it's
#    always designated (1 or 2).
#

# Step 1: Define truth values
F, G, T = 0, 1, 2
VALUES = {'F': F, 'G': G, 'T': T}
NAMES = {F: 'F', G: 'G', T: 'T'}

# Step 2: Define logical connectives
def op_not(a):
    if a == T:
        return F
    if a == F:
        return T
    return G  # NOT(G) = G

def op_and(a, b):
    return min(a, b)

def op_or(a, b):
    return max(a, b)

def is_designated(val):
    return val in [G, T]

# Step 3: Test the formula P = ¬(A ∧ B) ∨ (B ∧ A)
def check_formula_l():
    """
    Checks if the core formula from argument L is a tautology.
    The argument is A |- (A & B) -> (B & A), interpreted as A |- BOX(NOT(A&B) OR (B&A)).
    This is valid if P = NOT(A&B) OR (B&A) is a tautology.
    """
    print("Testing if the formula P = ¬(A ∧ B) ∨ (B ∧ A) is a tautology in 3-valued logic.")
    print("-" * 50)
    print("A | B | A ∧ B | B ∧ A | ¬(A ∧ B) | P = ¬(A ∧ B) ∨ (B ∧ A) | Is Designated?")
    print("-" * 50)

    all_tautological = True
    for val_a_name, val_a in VALUES.items():
        for val_b_name, val_b in VALUES.items():
            # Calculate sub-expressions
            a_and_b = op_and(val_a, val_b)
            b_and_a = op_and(val_b, val_a)
            not_a_and_b = op_not(a_and_b)
            
            # Final result for the formula P
            final_value = op_or(not_a_and_b, b_and_a)

            designated_status = is_designated(final_value)
            if not designated_status:
                all_tautological = False

            # Output the step-by-step calculation for this row
            print(
                f"{val_a_name} | {val_b_name} | "
                f"{NAMES[a_and_b]:^5} | {NAMES[b_and_a]:^5} | "
                f"{NAMES[not_a_and_b]:^8} | "
                f"{NAMES[final_value]:^26} | " # Show final result
                f"{designated_status}"
            )
            # Example of printing each number in the equation for one case
            if val_a_name == 'G' and val_b_name == 'T':
                print(f"  └─ Example Calc: ¬({val_a} ∧ {val_b}) ∨ ({val_b} ∧ {val_a})")
                print(f"                 = ¬({a_and_b}) ∨ ({b_and_a})")
                print(f"                 = {not_a_and_b} ∨ {b_and_a}")
                print(f"                 = {final_value}")

    print("-" * 50)
    if all_tautological:
        print("\nConclusion: The formula is a tautology.")
        print("Since the formula inside the box `□` is a tautology, the conclusion of argument L is a logical truth.")
        print("An argument with a logically true conclusion is always valid.")
    else:
        print("\nConclusion: The formula is NOT a tautology.")

if __name__ == "__main__":
    check_formula_l()
