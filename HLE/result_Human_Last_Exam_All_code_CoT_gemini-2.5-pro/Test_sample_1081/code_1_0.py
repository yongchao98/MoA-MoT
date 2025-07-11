def solve_e8_tori_number():
    """
    This script calculates and displays the formula for the number of
    F_q-rational maximal tori of a reductive group of type E_8.
    """
    # Parameters for the reductive group G of type E_8
    group_type = "E_8"
    rank = 8
    num_positive_roots = 120
    degrees = [2, 8, 12, 14, 18, 20, 24, 30]

    # Construct the formula string
    # Numerator part: |G(F_q)|
    product_terms = []
    for d in degrees:
        product_terms.append(f"(q^{d} - 1)")
    
    numerator_str = f"q^{num_positive_roots} * {' * '.join(product_terms)}"

    # Denominator part: |T_0(F_q)|
    denominator_str = f"(q - 1)^{rank}"

    # Print the final equation
    print(f"The exact number of F_q-rational maximal tori of a reductive group G of type {group_type} is given by the formula N, where q is the number of elements in the field F_q:")
    print(f"N = ({numerator_str}) / ({denominator_str})")
    print("\n" + "="*80 + "\n")

    # Output each number in the final equation as requested
    print("The specific numerical parameters in this formula for type E_8 are:")
    print(f"- The rank of the group: {rank}")
    print(f"- The number of positive roots (exponent of q in the leading factor): {num_positive_roots}")
    
    degree_list_str = ', '.join(map(str, degrees))
    print(f"- The degrees of the fundamental invariants (exponents in the polynomial factors): {degree_list_str}")

solve_e8_tori_number()