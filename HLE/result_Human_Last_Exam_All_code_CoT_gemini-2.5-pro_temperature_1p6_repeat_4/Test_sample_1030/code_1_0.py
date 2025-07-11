# Helper function to map numeric values to truth value names
def val_to_name(v):
    if v == 1:
        return 'T (1)'
    elif v == 0.5:
        return 'G (0.5)'
    return 'F (0)'

# Define the logical connectives for a 3-valued logic (LP)
def negation(v):
    return 1 - v

def conjunction(v1, v2):
    return min(v1, v2)
    
def disjunction(v1, v2):
    return max(v1, v2)

def implication(v1, v2):
    # Standard material implication: A -> B is defined as not(A) or B
    return disjunction(negation(v1), v2)

# Define the set of truth values and designated values
values = [1, 0.5, 0]  # T, G, F
designated = {1, 0.5}

# Test the conclusion of argument L: (A ∧ B) → (B ∧ A)
def check_l_conclusion():
    print("Evaluating conclusion of argument L: (A ∧ B) → (B ∧ A)")
    print("-" * 50)
    print(f"{'v(A)':<10} | {'v(B)':<10} | {'v(A ∧ B)':<12} | {'v(B ∧ A)':<12} | {'v(Conclusion)':<15} | {'Designated?':<12}")
    print("-" * 75)
    
    is_tautology = True
    for v_a in values:
        for v_b in values:
            # Evaluate sub-expressions
            val_a_and_b = conjunction(v_a, v_b)
            val_b_and_a = conjunction(v_b, v_a)
            
            # Evaluate the final conclusion
            final_val = implication(val_a_and_b, val_b_and_a)
            
            is_designated = "Yes" if final_val in designated else "No"
            if not (final_val in designated):
                is_tautology = False

            # Print the truth table row
            print(f"{val_to_name(v_a):<10} | {val_to_name(v_b):<10} | {val_to_name(val_a_and_b):<12} | {val_to_name(val_b_and_a):<12} | {val_to_name(final_val):<15} | {is_designated:<12}")
            
    print("-" * 75)
    if is_tautology:
        print("\nResult: The conclusion is a tautology (always designated).")
        print("Since the conclusion is always true, the argument 'A ⊢ (A ∧ B) → (B ∧ A)' is valid regardless of the premise.")
    else:
        print("\nResult: The conclusion is NOT a tautology.")

check_l_conclusion()