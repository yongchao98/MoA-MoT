import pygap

def solve():
    """
    Computes the Schur multiplier of a given permutation group G,
    and then counts the number of its proper subgroups up to isomorphism.
    """
    print("This script uses the 'pygap' library to interface with the GAP system.")
    print("Please ensure GAP is installed on your system for this script to run.")
    
    # Initialize GAP
    try:
        gap = pygap.gap_instance(max_workspace_size=2048)
    except Exception as e:
        print("\nFailed to initialize GAP. The calculation cannot proceed.")
        print(f"Error: {e}")
        return

    # The problem defines a group on the set {1, ..., 9, x, y, z}.
    # We map the symbols x, y, z to numbers 10, 11, 12 for use in GAP.
    a_str = "(1, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10)"  # y=11, x=10
    b_str = "(1, 8, 5, 9)(4, 10, 7, 6)"            # x=10
    c_str = "(1, 2)(3, 12)(4, 8)(5, 6)(7, 11)(9, 10)" # z=12, y=11, x=10
    
    print("\nStep 1: Defining the group G from its generators in GAP.")
    G = gap.eval(f"Group({a_str}, {b_str}, {c_str})")
    print("Group G has been created.")

    print("\nStep 2: Computing the Schur multiplier A = M(G).")
    schur_multiplier_A = gap.SchurMultiplier(G)
    invariants = list(map(int, schur_multiplier_A.AbelianInvariants()))

    # Handle the case of a trivial multiplier
    if not invariants or (len(invariants) == 1 and invariants[0] == 1):
        invariants = [] # Canonical representation for trivial group
        A_structure = "C_1 (the trivial group)"
        invariants_tuple = tuple([])
    else:
        A_structure_parts = [f"Z_{n}" for n in invariants]
        A_structure = " x ".join(A_structure_parts)
    
    print(f"The Schur multiplier A is isomorphic to the abelian group {A_structure}.")

    # To find the number of subgroup isomorphism classes, we create a model group for A
    # in GAP and ask it to find all subgroups and their structures.
    print("\nStep 3: Finding all subgroup isomorphism types for A.")
    A_model = gap.AbelianGroup(invariants)
    subgroups = A_model.Subgroups()
    
    isomorphism_types = set()
    for sub in subgroups:
        # The isomorphism type is determined by the abelian invariants of the subgroup.
        sub_invariants = tuple(sorted(list(map(int, sub.AbelianInvariants()))))
        isomorphism_types.add(sub_invariants)
        
    total_iso_classes = len(isomorphism_types)
    num_proper_iso_classes = total_iso_classes - 1
    
    # Present the found types for clarity
    type_strs = []
    for t in sorted(list(isomorphism_types), key=lambda x: (len(x), x)):
        if not t:
            type_strs.append("C_1")
        else:
            type_strs.append(" x ".join([f"Z_{i}" for i in t]))

    print(f"The group A has a total of {total_iso_classes} distinct subgroup isomorphism types.")
    print(f"The found types are: {type_strs}")

    print("\nA proper subgroup is any subgroup except the group A itself.")
    print("The number of isomorphism classes of proper subgroups is found by subtracting 1 from the total count.")

    print("\nFinal Calculation:")
    print(f"Total isomorphism types - 1 = Number of proper isomorphism types")
    print(f"{total_iso_classes} - 1 = {num_proper_iso_classes}")
    
    print(f"\nThe number of proper subgroups of A, up to isomorphism, is {num_proper_iso_classes}.")

if __name__ == '__main__':
    solve()