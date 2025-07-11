def get_number_of_tori_e8():
    """
    Calculates the number of rational maximal tori of a group of type E8
    over a finite field Fq.
    """
    # The total number of roots in the E8 root system is 240.
    total_roots_e8 = 240

    # The number of positive roots, N, is half of the total number of roots.
    num_positive_roots = total_roots_e8 // 2

    # The number of Fq-rational maximal tori is given by the formula q^N.
    # The problem asks for the final equation and its numerical components.
    
    print("The exact number of rational maximal tori for a group of type E8 over Fq is given by the equation:")
    # We build the equation string. The base is the variable 'q'. The exponent is the number we calculated.
    equation_str = f"Number = q^{num_positive_roots}"
    
    print(equation_str)
    
    # We output the single numerical value that appears in the equation as requested.
    print(f"\nThe only number in this equation is the exponent: {num_positive_roots}")

get_number_of_tori_e8()