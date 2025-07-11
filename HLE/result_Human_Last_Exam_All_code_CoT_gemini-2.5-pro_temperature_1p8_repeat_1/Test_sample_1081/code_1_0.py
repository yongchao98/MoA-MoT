def solve_tori_count():
    """
    Calculates the number of F_q-rational maximal tori of a reductive group of type E_8.

    The number of F_q-rational maximal tori of a reductive group G over F_q
    is given by the number of F-conjugacy classes of its Weyl group W.
    
    For a group G of type E_8, the Dynkin diagram has trivial automorphism group.
    This implies that any such group G over F_q is split.

    For a split group, the Frobenius endomorphism acts trivially on the Weyl group.
    Thus, F-conjugacy classes are the same as ordinary conjugacy classes.

    The problem reduces to finding the number of conjugacy classes of the Weyl group W(E_8).
    This is a known result from the theory of finite reflection groups.
    """

    # Number of conjugacy classes in the Weyl group W(E_8)
    num_classes_W_E8 = 112

    print("Let N be the number of F_q-rational maximal tori of a reductive group of type E_8.")
    print("N is equal to the number of conjugacy classes of the corresponding Weyl group, W(E_8).")
    
    # The final calculation is just stating this known value.
    print("\nThe number of conjugacy classes of W(E_8) is:")
    print(num_classes_W_E8)
    
    print("\nTherefore, the final equation for the number of tori is:")
    print(f"N = {num_classes_W_E8}")

solve_tori_count()