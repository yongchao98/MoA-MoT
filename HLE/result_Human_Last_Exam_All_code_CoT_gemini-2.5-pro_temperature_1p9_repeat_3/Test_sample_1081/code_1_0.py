def calculate_e8_tori_count():
    """
    This function calculates the number of F_q-rational maximal tori
    of a reductive group of type E_8 over a finite field F_q.

    The number of such tori is equal to the order of the Weyl group W(E_8).
    This function computes |W(E_8)| from its prime factorization.
    """

    # The prime factorization of the order of the Weyl group of type E_8
    # is |W(E_8)| = 2^14 * 3^5 * 5^2 * 7^1.
    factors = {
        2: 14,
        3: 5,
        5: 2,
        7: 1
    }

    # Calculate each term in the product
    term_values = {base: base**exp for base, exp in factors.items()}
    
    # Calculate the total order
    total_order = 1
    for val in term_values.values():
        total_order *= val

    # Print the explanation and the step-by-step calculation
    print("The number of F_q-rational maximal tori of a group of type E_8 is given by the order of its Weyl group, |W(E_8)|.")
    print("The order is calculated from its known prime factorization.\n")
    
    equation_str = " * ".join([f"{base}^{exp}" for base, exp in factors.items()])
    print(f"|W(E_8)| = {equation_str}\n")
    
    print("Calculating each term:")
    value_str_list = []
    for base, exp in factors.items():
        value = term_values[base]
        print(f"{base}^{exp} = {value}")
        value_str_list.append(str(value))
    
    print("\nPutting it all together:")
    final_equation_str = " * ".join(value_str_list)
    print(f"|W(E_8)| = {final_equation_str} = {total_order}")

# Execute the function to print the result
calculate_e8_tori_count()