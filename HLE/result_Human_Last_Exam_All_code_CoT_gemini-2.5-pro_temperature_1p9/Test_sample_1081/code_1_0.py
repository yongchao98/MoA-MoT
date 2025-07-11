def solve():
    """
    This function calculates the number of F_q-rational maximal tori of a
    reductive group G of type E_8 over the finite field F_q.
    
    The number of F_q-rational maximal tori (up to F_q-conjugacy) in a split
    reductive group is equal to the number of conjugacy classes of its Weyl group.
    For a group of type E_8, the Weyl group is W(E_8).
    
    The number of conjugacy classes of the Weyl group W(E_8) is a known
    mathematical fact, which is 112.
    """
    
    # The number of conjugacy classes of the Weyl group W(E_8).
    num_tori = 112
    
    # Print the explanation and the final equation.
    print("The number of F_q-rational maximal tori is equal to the number of conjugacy classes of the Weyl group W(E_8).")
    print(f"Number of rational maximal tori = Number of conjugacy classes of W(E_8) = {num_tori}")

solve()