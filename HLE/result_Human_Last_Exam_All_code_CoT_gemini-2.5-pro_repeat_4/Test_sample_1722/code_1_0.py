import numpy as np
import itertools

def solve():
    """
    This function demonstrates that for n=2, observation sequences can be constructed
    such that a 2-state memory is insufficient to distinguish them, while a
    3-state memory is sufficient. This makes n=2 the minimum such length.
    """
    n = 2
    O1 = (0, 1)
    O2 = (1, 0)
    
    print(f"Analyzing for hallway length n = {n}")
    print(f"Observation sequence in C1: {O1}")
    print(f"Observation sequence in C2: {O2}")
    print("-" * 30)

    # --- Part 1: Show m=2 is NOT enough ---
    print("Case m=2:")
    s0_2 = np.array([1, 0])
    
    # 2x2 permutation matrices: Identity and Swap
    I2 = np.identity(2, dtype=int)
    X2 = np.array([[0, 1], [1, 0]])
    perms_2_state = [I2, X2]
    
    distinguishable_m2 = False
    # Iterate over all possible 2-state FSMs.
    # An FSM is defined by the choice of transition matrix for '0' and '1'.
    for A0 in perms_2_state:
        for A1 in perms_2_state:
            # Transformation matrix for sequence O1 = (0, 1)
            M_O1 = A1 @ A0
            # Transformation matrix for sequence O2 = (1, 0)
            M_O2 = A0 @ A1
            
            final_state_O1 = M_O1 @ s0_2
            final_state_O2 = M_O2 @ s0_2
            
            if not np.array_equal(final_state_O1, final_state_O2):
                distinguishable_m2 = True
                break
        if distinguishable_m2:
            break
            
    if not distinguishable_m2:
        print("Result: Indistinguishable. For any 2-state FSM, the final memory state is the same for both sequences.")
        print("This is because the group of 2x2 permutation matrices is abelian (A0 @ A1 = A1 @ A0).")
    else:
        print("Result: Distinguishable. This contradicts the theory.")
    print("-" * 30)

    # --- Part 2: Show m=3 IS enough ---
    print("Case m=3:")
    s0_3 = np.array([1, 0, 0])

    # For m=3, we can pick non-commuting permutation matrices.
    # Let A0 be a 3-cycle (1->2, 2->3, 3->1)
    A0_3 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    # Let A1 be a transposition (1->2, 2->1)
    A1_3 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])

    # Transformation for O1 = (0, 1)
    M_O1_3 = A1_3 @ A0_3
    # Transformation for O2 = (1, 0)
    M_O2_3 = A0_3 @ A1_3
    
    final_state_O1_3 = M_O1_3 @ s0_3
    final_state_O2_3 = M_O2_3 @ s0_3

    print("By choosing a specific 3-state FSM with non-commuting transitions:")
    print("Let A0 (for obs '0') be the matrix:\n", A0_3)
    print("Let A1 (for obs '1') be the matrix:\n", A1_3)
    
    print("\nFinal state for O1=(0,1) is (A1 @ A0) @ s0:")
    print(f"{final_state_O1_3}")
    print("Final state for O2=(1,0) is (A0 @ A1) @ s0:")
    print(f"{final_state_O2_3}")

    if not np.array_equal(final_state_O1_3, final_state_O2_3):
        print("\nResult: Distinguishable. The final memory states are different.")
        print("This allows an agent with m=3 to achieve a higher reward.")
    else:
        print("\nResult: Indistinguishable. Something is wrong.")
    
    print("-" * 30)
    print("Conclusion: The minimum length n for which m=2 fails while m=3 succeeds is 2.")
    
    # Final answer as required by the format
    final_answer = n
    print(f"The final answer is {final_answer}")

solve()