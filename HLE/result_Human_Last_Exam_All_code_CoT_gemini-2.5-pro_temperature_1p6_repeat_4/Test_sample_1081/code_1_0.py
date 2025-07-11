def solve_e8_tori_problem():
    """
    Calculates the number of Fq-rational maximal tori of a reductive group of type E8.
    
    The problem is solved based on a key theorem in the theory of reductive groups.
    """

    # Step 1: Theoretical connection
    # The number of Fq-rational maximal tori of a split reductive group G over a finite field Fq
    # is equal to the number of conjugacy classes of its Weyl group W.
    # For a group G of type E8, the Weyl group is W(E8).
    
    # Let N be the number of Fq-rational maximal tori.
    # The "equation" is: N = Number of conjugacy classes of W(E8)
    
    # Step 2: Look up the known mathematical value
    # The number of conjugacy classes of the Weyl group W(E8) is a well-established result
    # in mathematics, found in sources like R.W. Carter's paper "Conjugacy classes in the Weyl group".
    
    num_conjugacy_classes_W_E8 = 112
    
    # Step 3: Print the result and the "equation"
    # We output the components of the conceptual equation leading to the answer.
    equation_lhs = "Number of F_q-rational maximal tori of G(E_8)"
    equation_mid = "Number of conjugacy classes of W(E_8)"
    final_value = num_conjugacy_classes_W_E8

    print(f"The number of F_q-rational maximal tori is determined by the structure of the Weyl group.")
    print(f"The equation is: {equation_lhs} = {equation_mid}")
    print(f"The number of conjugacy classes for the Weyl group of type E_8 is a known value.")
    print(f"The number in the final equation is: {final_value}")
    print(f"Therefore, the total number of F_q-rational maximal tori is {final_value}.")

solve_e8_tori_problem()