def illustrate_contradiction():
    """
    This code illustrates the contradiction at the heart of the argument.
    The argument shows that if such a bounded uncountable set of functions X existed,
    it would imply the existence of a strictly increasing sequence of ordinals
    of length omega_1, with all its elements bounded by a single ordinal
    from omega_1. This is impossible.

    We use non-negative integers as an analogy for countable ordinals.
    A sequence of length omega_1 is analogous to an infinite sequence of integers.
    The supposed bound `g(gamma)` is represented by a fixed integer `B`.
    """

    B = 100 # Our analogous g(gamma), a fixed upper bound.

    print(f"Let's assume there is an infinite, strictly increasing sequence of non-negative integers, `a_i`, bounded by B={B}.")
    print("The sequence is `a_0, a_1, a_2, ...` where `a_{i+1} > a_i`.")
    print(f"The assumption is that `a_i < B` for all `i`.")
    print("\nLet's prove by induction that `a_i >= i` for all `i`.")
    
    print("Base Case (i=0): `a_0` is a non-negative integer, so `a_0 >= 0`. This holds.")
    
    print("Inductive Step: Assume `a_k >= k` for some integer `k >= 0`.")
    print("Since the sequence is strictly increasing, `a_{k+1} > a_k`.")
    print("Because `a_k` and `a_{k+1}` are integers, `a_{k+1} >= a_k + 1`.")
    print(f"From our assumption `a_k >= k`, we have `a_k + 1 >= k + 1`.")
    print(f"Combining these, we get `a_{k+1} >= k + 1`.")
    print("The induction is complete. We have shown `a_i >= i` for all `i`.")

    print(f"\nNow, let's consider the element `a_B` of the sequence (i.e., the element at index {B}).")
    i = B
    print(f"According to our proof, `a_i >= i`, which means `a_{B} >= {B}`.")

    print(f"However, our initial assumption was that *all* elements are less than {B}, which implies `a_{B} < {B}`.")
    
    print("\nThis leads to a contradiction:")
    # The key contradiction expressed with the values.
    # The smallest possible value for a_B according to the induction.
    a_B_val = B
    
    # Printing each number in the final equation.
    print(f"The bound B = {B}.")
    print(f"The sequence element at index B is a_B.")
    print(f"From `a_i >= i`, we have a_B >= {B}.")
    print(f"From the boundedness assumption, we have a_B < {B}.")
    
    print("\nFinal contradiction shown with values:")
    print(f"{a_B_val} >= {B}  AND  {a_B_val} < {B}")
    
    print("\nThis is impossible. The initial assumption must be false.")
    print("Therefore, an infinite strictly increasing sequence of integers cannot be bounded.")
    print("Analogously, a strictly increasing sequence of ordinals of length omega_1 cannot be bounded by an ordinal in omega_1.")


illustrate_contradiction()