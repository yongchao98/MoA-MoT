def explain_cohomology_degrees():
    """
    This function explains the significance of low-degree cohomology groups
    in semi-abelian categories to answer the user's question.
    """
    print("Step-by-step analysis of cohomology degrees:")
    print("============================================")

    # Degree 0
    print("Degree 0 (H^0):")
    print("The 0-th cohomology group, H^0(B, M), represents the B-invariants of the module M.")
    print("This degree deals with fixed points, not extensions or obstructions.")
    print("-" * 40)

    # Degree 1
    print("Degree 1 (H^1):")
    print("The 1st cohomology group, H^1(B, M), classifies derivations and is related to Ext^1(B, M), which classifies extensions of modules.")
    print("This is the minimal degree where 'non-trivial extensions' become significant.")
    print("However, the 'obstruction' concept is more clearly and fundamentally expressed in the next degree.")
    print("-" * 40)
    
    # Degree 2
    print("Degree 2 (H^2):")
    print("The 2nd cohomology group, H^2(B, M), is pivotal. It classifies extensions of the object B by the module M (a generalization of group extensions).")
    print("A non-zero element in H^2(B, M) serves two roles simultaneously:")
    print("  1. It represents a class of 'non-trivial extensions'.")
    print("  2. It acts as the 'obstruction' to the extension being split (i.e., being equivalent to a semi-direct product).")
    print("Therefore, degree 2 is the minimal degree where the concepts of extensions AND obstructions are both fundamentally significant and intertwined.")
    print("-" * 40)

    # Degree 3
    print("Degree 3 (H^3):")
    print("The 3rd cohomology group, H^3(B, M), is typically associated with obstructions to extending more complex structures, such as crossed modules.")
    print("While it is about obstructions, it is not the *minimal* degree where the interplay with extensions begins.")
    print("-" * 40)
    
    print("\nConclusion:")
    print("The minimal cohomology degree where both non-trivial extensions and obstructions are significant is 2.")

# Execute the explanation
explain_cohomology_degrees()

final_answer = "C"
print(f"\nThe corresponding answer choice for degree 2 is {final_answer}.")