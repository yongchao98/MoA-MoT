def solve_e8_tori_count():
    """
    Calculates the number of F_q-rational maximal tori of a split reductive group of type E_8.

    The number of F_q-rational maximal tori in a split reductive group G is equal to
    the number of conjugacy classes in its Weyl group W. For a group of type E_8, the
    Weyl group is W(E_8). The number of conjugacy classes of W(E_8) is a known
    mathematical constant.
    """
    
    # The number of conjugacy classes of the Weyl group of type E8.
    num_conjugacy_classes_W_E8 = 112
    
    print("Let G be a split reductive group of type E_8 over the finite field F_q.")
    print("The number of F_q-rational maximal tori of G is equal to the number of conjugacy classes of its Weyl group, W(E_8).")
    print("The number of conjugacy classes of W(E_8) is known to be:")
    print(num_conjugacy_classes_W_E8)

solve_e8_tori_count()