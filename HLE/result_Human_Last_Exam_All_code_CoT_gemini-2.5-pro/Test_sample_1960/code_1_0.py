def solve_equipartition_logic():
    """
    This function defines and prints the linear logic formulas for the equipartitioning problem.
    Since the formulas themselves are the answer, we will print their structure.
    """

    # Define the base formulas A and B without literals
    # B = (1 ⊸ ⊥)
    # A = (B ⊸ B) which is ((1 ⊸ ⊥) ⊸ (1 ⊸ ⊥))
    B_str = "(1 ⊸ ⊥)"
    A_str = f"({B_str} ⊸ {B_str})"

    print("The solution uses two base formulas, A and B, constructed without literals:")
    print(f"Let B = {B_str}")
    print(f"Let A = {A_str}\n")

    print("The function f(w) is defined as:")
    # f(w) = (A ⊸ B^w) ⊗ A
    # We use 'B^w' as a shorthand for B ⊗ B ⊗ ... ⊗ B (w times)
    print("f(w) = (A ⊸ (B ⊗ B ⊗ ... ⊗ B (w times))) ⊗ A")
    print(f"f(w) = ({A_str} ⊸ ({B_str}^w)) ⊗ {A_str}\n")
    
    print("The goal formula C(W, m, b) is defined as:")
    # C = (B^b)^m
    print("C = (B ⊗ ... ⊗ B (b times)) ⊗ ... ⊗ (B ⊗ ... ⊗ B (b times))  (m groups)")
    print(f"C = ({B_str}^b)^m\n")

    print("The sequent to be proved is:")
    print("{ f(w) | w ∈ W } vdash C")
    print("Which translates to:")
    print(f"{{ (({A_str} ⊸ ({B_str}^w)) ⊗ {A_str}) | w ∈ W }} vdash (({B_str}^b)^m)")


solve_equipartition_logic()
