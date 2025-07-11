import numpy as np

def get_extension_group_table(cohomology_class):
    """
    Constructs the multiplication table for a group extension of Z_2 by Z_2.
    The structure is determined by a class in H^2(Z_2, Z_2).

    Args:
        cohomology_class (int): 0 for the trivial class, 1 for the non-trivial class.
                                This value comes from the 2nd cohomology group.
    """
    B = [0, 1]  # The quotient group Z_2
    A = [0, 1]  # The kernel group Z_2
    
    # The group elements are pairs (a, b) with a in A, b in B.
    # We map them to integers 0, 1, 2, 3 for easier indexing.
    # (0,0)->0, (1,0)->1, (0,1)->2, (1,1)->3
    elements = [(a, b) for b in B for a in A]
    
    # The multiplication is (a1, b1) * (a2, b2) = (a1 + a2 + f(b1, b2), b1 + b2)
    # where f is the 2-cocycle. All additions are mod 2.
    
    # Define the 2-cocycle f(b1, b2) based on the cohomology class.
    # The action of B on A is trivial.
    if cohomology_class == 0:
        # Trivial cocycle f(b1, b2) = 0. This corresponds to the split extension.
        # The resulting group is Z_2 x Z_2.
        f = lambda b1, b2: 0
        group_name = "Z_2 x Z_2 (trivial extension)"
    elif cohomology_class == 1:
        # Non-trivial cocycle f(1, 1) = 1, and 0 otherwise.
        # This corresponds to the non-split extension.
        # The resulting group is Z_4.
        f = lambda b1, b2: 1 if b1 == 1 and b2 == 1 else 0
        group_name = "Z_4 (non-trivial extension)"
    else:
        print("Invalid cohomology class. Please use 0 or 1.")
        return

    table = np.zeros((4, 4), dtype=int)
    for i, (a1, b1) in enumerate(elements):
        for j, (a2, b2) in enumerate(elements):
            # Apply the group law
            b_res = (b1 + b2) % 2
            a_res = (a1 + a2 + f(b1, b2)) % 2
            
            # Find the index of the resulting element
            res_element = (a_res, b_res)
            res_index = elements.index(res_element)
            table[i, j] = res_index

    print(f"This calculation demonstrates how a class from the 2nd cohomology group determines the extension structure.")
    print(f"The minimal degree at which non-trivial extensions become significant is 2.")
    print(f"Final equation for the element operation: (a_res, b_res) = ((a1 + a2 + f(b1, b2)) mod 2, (b1 + b2) mod 2)")
    print(f"\nUsing cohomology_class = {cohomology_class}, the resulting group is {group_name}.")
    print("Group elements map to integers: (0,0)->0, (1,0)->1, (0,1)->2, (1,1)->3")
    print("Multiplication Table:")
    print(table)


# Run for the trivial class (split extension)
get_extension_group_table(0)
print("\n" + "="*50 + "\n")
# Run for the non-trivial class (non-split extension)
get_extension_group_table(1)
