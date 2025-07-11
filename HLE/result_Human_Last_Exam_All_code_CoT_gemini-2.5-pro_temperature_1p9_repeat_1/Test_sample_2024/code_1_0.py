import sys

# Suppress warnings if any might occur with Unicode printing in certain environments
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

def explain_logic():
    """Prints the step-by-step logical deduction."""
    print("### Step-by-Step Derivation ###")
    print("\n1. Analyze the Axiom Truth Value and its consequences:")
    print("   Axiom: ∀x∀y∀z (T(x,y,z) → □(∀w (R(z,w) → T(x,y,w))))")
    print("   Given: The accessibility relation R is an equivalence relation on {w1, w2, w3},")
    print("   meaning all worlds are mutually accessible.")
    print("   Given: The system has truth values {0, 0.5, 1} and 'T(0.5, a, w1) holds'.")
    print("\n   Consequence (a): The value of any predicate T(x,y,z) cannot be 0.5.")
    print("   - If val(T(x,y,z)) were 0.5, the axiom forces its consequent to have value 1.")
    print("   - This consequent, due to R, implies val(T(x,y,w')) = 1 for ALL worlds w'.")
    print("   - This includes the original T(x,y,z), meaning its value must be 1. This contradicts the assumption that its value was 0.5.")
    print("\n   Consequence (b): The value of T(x,y,z) is rigid (constant) across all worlds.")
    print("   - If val(T(x,y,z), w_k) = 1 for any world w_k, the axiom forces val(T(x,y,z), w_j) = 1 for ALL worlds w_j.")
    print("\n   Therefore, for any instance (x,y,z), its truth value, v, is in {0, 1} and is constant across all worlds.")

def evaluate_statement():
    """Evaluates the statement bottom-up and prints the values."""
    print("\n2. Evaluate the statement bottom-up: □(∀x∀y∀z (T(x,y,z) → □(T(x,y,z))))")
    
    # We analyze the value of the implication I = T -> []T for any instance.
    # Let v be the value of T(x,y,z), which can be 0 or 1.
    
    # Case v = 0:
    val_T_case_0 = 0
    val_inner_Box_case_0 = 0 # min(0, 0, 0)
    val_implication_case_0 = max(1 - val_T_case_0, val_inner_Box_case_0) # max(1, 0) = 1
    
    # Case v = 1:
    val_T_case_1 = 1
    val_inner_Box_case_1 = 1 # min(1, 1, 1)
    val_implication_case_1 = max(1 - val_T_case_1, val_inner_Box_case_1) # max(0, 1) = 1
    
    # Since the implication is 1 in all cases, the universal quantifier results in 1.
    val_forall = 1
    
    # The outer box operator on a constant value of 1 results in 1.
    val_outer_box = 1

    print("\n   - Let v be the value of T(x,y,z). v must be 0 or 1.")
    print(f"   - The value of the inner □(T(x,y,z)) is min(v, v, v) = v.")
    print(f"   - The value of the implication T(x,y,z) → □(T(x,y,z)) is max(1-v, v).")
    print(f"     - If v = 0, the value is max(1-0, 0) = {val_implication_case_0}.")
    print(f"     - If v = 1, the value is max(1-1, 1) = {val_implication_case_1}.")
    print("     Thus, the implication always has a value of 1.")
    
    print("\n   - The value of ∀x∀y∀z (...) is the minimum over all instances.")
    print(f"     Since every instance yields 1, the result is {val_forall}.")

    print("\n   - The value of the outer □(...) is the minimum over all worlds.")
    print(f"     Since the inner statement has a constant value of 1, the result is {val_outer_box}.")

    return val_outer_box

def print_final_equation_values():
    """Prints the final equation with the determined values for each part."""
    print("\n### Final Evaluation ###")
    print("Expression: □ ( ∀x∀y∀z ( T(x,y,z) → □(T(x,y,z)) ) )")
    print("Values:     1 (    1      (     v     → 1(    v   ) ) )")
    print("                     where v = {0, 1}, and for both values, (v → v) = 1")
    
    final_value_inner_T = "{0, 1}"
    final_value_inner_box = "v"
    final_value_implication = "1"
    final_value_forall = "1"
    final_value_outer_box = "1"

    print("\nBreaking down the equation with final values:")
    print(f"1. T(x, y, z): has a consistent value v, which is either 0 or {val_implication_case_1}.")
    print(f"2. □(T(x, y, z)): The value is min(v,v,v), which is also v.")
    print(f"3. T(x, y, z) → □(T(x, y, z)): The value is max(1-v, v), which is always {final_value_implication}.")
    print(f"4. ∀x∀y∀z (...): The minimum over all instances, each of which is 1, results in {final_value_forall}.")
    print(f"5. □(...): The value is the minimum over all worlds for a statement of value 1, resulting in {final_value_outer_box}.")


# Main execution
explain_logic()
final_answer = evaluate_statement()
print_final_equation_values()
print(f"\nThe determined truth value of the statement is {final_answer}.")
print("<<<1>>>")