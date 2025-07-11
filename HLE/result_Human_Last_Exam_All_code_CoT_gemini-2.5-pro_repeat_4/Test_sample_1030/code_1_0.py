# Define truth values using integers for easy comparison
T, G, F = 2, 1, 0
val_map = {2: 'T', 1: 'G', 0: 'F'}
DESIGNATED_VALUES = {T, G}

# Define the 3-valued logic operations
def neg(v):
    if v == T: return F
    if v == F: return T
    return G  # v == G

def conj(v1, v2):
    return min(v1, v2)

def disj(v1, v2):
    return max(v1, v2)

def impl(v1, v2):
    return disj(neg(v1), v2)

def check_validity_K():
    """
    Checks the validity of the argument K: A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B)
    """
    print("Checking validity of Argument K: A ∧ B vdash (¬A ∨ ¬B) → (A ∧ B)")
    print("Truth Values: T=2, G=1, F=0")
    print("Designated Values: {T, G} or {2, 1}\n")

    is_valid = True
    # Iterate through all 3*3 = 9 possible truth value assignments for A and B
    for v_a in [T, G, F]:
        for v_b in [T, G, F]:
            # 1. Calculate the value of the premise: P = A ∧ B
            premise_val = conj(v_a, v_b)

            # 2. Calculate the value of the conclusion: C = (¬A ∨ ¬B) → (A ∧ B)
            # Breaking it down step-by-step for clarity
            neg_a = neg(v_a)
            neg_b = neg(v_b)
            left_of_arrow = disj(neg_a, neg_b)
            right_of_arrow = conj(v_a, v_b) # same as premise
            conclusion_val = impl(left_of_arrow, right_of_arrow)

            # 3. Check the validity condition
            # If premise is designated, conclusion must also be designated
            validity_check = "Holds"
            if premise_val in DESIGNATED_VALUES:
                if conclusion_val not in DESIGNATED_VALUES:
                    is_valid = False
                    validity_check = "Fails"
            else:
                validity_check = "Holds (Premise not designated)"


            # Print the detailed evaluation for this assignment
            print(f"Case: v(A)={val_map[v_a]}, v(B)={val_map[v_b]}")
            print(f"  - Premise: v(A ∧ B) = v({val_map[v_a]} ∧ {val_map[v_b]}) = {val_map[premise_val]}")
            print(f"  - Conclusion Part 1 (Left of →): v(¬A ∨ ¬B) = v(¬{val_map[v_a]} ∨ ¬{val_map[v_b]}) = v({val_map[neg_a]} ∨ {val_map[neg_b]}) = {val_map[left_of_arrow]}")
            print(f"  - Conclusion Part 2 (Right of →): v(A ∧ B) = {val_map[right_of_arrow]}")
            print(f"  - Conclusion (Full): v(({val_map[left_of_arrow]}) → ({val_map[right_of_arrow]})) = {val_map[conclusion_val]}")
            print(f"  - Validity Check: Premise is {val_map[premise_val]}, Conclusion is {val_map[conclusion_val]}. Rule: if Premise is T/G, Conclusion must be T/G. -> {validity_check}\n")


    if is_valid:
        print("="*40)
        print("Result: The argument K is VALID.")
        print("In all cases where the premise (A ∧ B) has a designated value (T or G),")
        print("the conclusion ((¬A ∨ ¬B) → (A ∧ B)) also has a designated value.")
        print("="*40)
    else:
        print("="*40)
        print("Result: The argument K is INVALID.")
        print("="*40)

if __name__ == '__main__':
    check_validity_K()
