import itertools

def check_non_group_of_size_2():
    """
    This function demonstrates that a set of size 2 can form a non-group structure.
    It checks the associativity of the logical NAND operation on the set {True, False}.
    """
    # Let our set G have two elements. We use {True, False} for convenience.
    G = [True, False]
    n = len(G)
    print(f"Checking for a set G of size n = {n}. G = {{True, False}}")
    print("The binary operation is logical NAND: a . b = NOT (a AND b)\n")

    # Define the NAND operation
    def nand_op(a, b):
        return not (a and b)

    print("Checking the associativity axiom: (a . b) . c = a . (b . c)")
    
    found_counterexample = False
    # Iterate through all possible combinations of a, b, c from G
    for a, b, c in itertools.product(G, repeat=3):
        left_side = nand_op(nand_op(a, b), c)
        right_side = nand_op(a, nand_op(b, c))

        if left_side != right_side:
            print(f"Found a counterexample for a={a}, b={b}, c={c}:")
            
            # Show the calculation step-by-step
            # Left side
            step1_left = nand_op(a, b)
            print(f"  Left side: ({a} NAND {b}) NAND {c}")
            print(f"           = {step1_left} NAND {c}")
            print(f"           = {left_side}")

            # Right side
            step1_right = nand_op(b, c)
            print(f"  Right side: {a} NAND ({b} NAND {c})")
            print(f"            = {a} NAND {step1_right}")
            print(f"            = {right_side}")
            
            print(f"\nSince {left_side} != {right_side}, the operation is not associative.")
            found_counterexample = True
            break
    
    print("\n--- Conclusion ---")
    print("Any set of size n=1 with a binary operation forms a group (the trivial group).")
    print(f"We have demonstrated that a set of size n=2 can have a binary operation that is not a group.")
    print("Therefore, the smallest number n (for a non-empty set) is 2.")


check_non_group_of_size_2()
