def print_upper_bound_formula():
    """
    This function prints the general formula for the upper bound on the maximum norm
    (k_{k,inf}) for the lattice of integers in a number field.
    """
    
    # The numbers that appear in the formula are 4 and 1.
    # As requested, we will output these numbers in the final equation.
    number_four = 4
    number_one = 1
    
    # The formula provides an upper bound for k_{k,inf} in terms of:
    # V: the covolume of the lattice
    # n: the degree of the number field
    # s: the number of pairs of complex embeddings of the field
    # pi: the mathematical constant pi
    
    # The inequality is: k_{k,inf} <= ( (4/pi)^s * V )^(1/n)
    
    print("The general upper bound for the maximum norm (k_{k,inf}) in relation to the covolume (V) is given by the following inequality:")
    print(f"k_{{k,inf}} <= (({number_four}/pi)^s * V)^({number_one}/n)")

print_upper_bound_formula()