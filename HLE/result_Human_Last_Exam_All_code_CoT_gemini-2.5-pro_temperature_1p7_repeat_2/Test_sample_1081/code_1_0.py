def solve_tori_question():
    """
    Calculates the number of F_q-rational maximal tori of a reductive group G of type E_8.

    The number of F_q-rational maximal tori of a split reductive group G over F_q
    is given by the number of conjugacy classes in its corresponding Weyl group W.
    """
    
    # The group G is of type E_8. The standard simple group of this type is split
    # over any field F_q.
    group_type = "E_8"
    
    # The corresponding Weyl group is W(E_8).
    weyl_group_type = "W(E_8)"
    
    # The problem reduces to finding the number of conjugacy classes of W(E_8).
    # This is a known mathematical fact.
    num_conjugacy_classes_W_E8 = 112
    
    # The "equation" is essentially an equality based on a theorem.
    # Number of tori = Number of conjugacy classes of W(E_8)
    
    print(f"Step 1: Identify the type of the reductive group G.")
    print(f"The group G is of type {group_type}.")
    print("-" * 30)
    
    print("Step 2: Relate the number of maximal tori to a property of the Weyl group.")
    print("For a split group G over a finite field, the number of F_q-rational maximal tori is equal")
    print(f"to the number of conjugacy classes of its Weyl group, which is {weyl_group_type}.")
    print("-" * 30)
    
    print(f"Step 3: Find the number of conjugacy classes for {weyl_group_type}.")
    print(f"The number of conjugacy classes of the Weyl group {weyl_group_type} is a known result in mathematics.")
    print(f"The final number in our 'equation' is: {num_conjugacy_classes_W_E8}")
    print("-" * 30)

    print("Conclusion:")
    print(f"The exact number of F_q-rational maximal tori of G is {num_conjugacy_classes_W_E8}.")

solve_tori_question()