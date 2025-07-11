import sys

def solve():
    """
    This function calculates the number of multiplication operations in the fully expanded
    expression for the sum bit s2 of a 3-bit binary adder.
    """
    # Step 1: Define the minimal sum-of-products for the carry bit c2.
    # c2 = (a1*b1) + (a1*c1) + (b1*c1), where c1 = a0*b0.
    # So, c2 = (a1*b1) + (a1*a0*b0) + (b1*a0*b0).
    # Each inner list represents a product term with its literals.
    c2_terms = [
        ["a1", "b1"],
        ["a1", "a0", "b0"],
        ["b1", "a0", "b0"]
    ]

    # Step 2: Define the minimal sum-of-products for the complement c2'.
    # c2' = a1'b1' + a1'c1' + b1'c1' where c1' = a0' + b0'.
    # So, c2' = a1'b1' + a1'(a0'+b0') + b1'(a0'+b0')
    # Expanding this gives: a1'b1' + a1'a0' + a1'b0' + b1'a0' + b1'b0'.
    c2_prime_terms = [
        ["a1'", "b1'"],
        ["a1'", "a0'"],
        ["a1'", "b0'"],
        ["b1'", "a0'"],
        ["b1'", "b0'"]
    ]

    # Step 3: Define the components of s2 based on its XOR expansion.
    # s2 = (a2'*b2'*c2) + (a2'*b2*c2') + (a2*b2'*c2') + (a2*b2*c2)
    s2_components = [
        ({"prefix": ["a2'", "b2'"], "name": "a2'b2'c2"}, c2_terms),
        ({"prefix": ["a2'", "b2"], "name": "a2'b2c2'"}, c2_prime_terms),
        ({"prefix": ["a2", "b2'"], "name": "a2b2'c2'"}, c2_prime_terms),
        ({"prefix": ["a2", "b2"], "name": "a2b2c2"}, c2_terms)
    ]

    total_multiplications = 0
    multiplication_counts = []
    
    print("The fully expanded expression for s2 is derived by substituting c2 and c2' into:")
    print("s2 = (a2' * b2' * c2) + (a2' * b2 * c2') + (a2 * b2' * c2') + (a2 * b2 * c2)\n")

    # Step 4 & 5: Expand each component and count multiplications.
    final_s2_expression_parts = []
    for component_info, term_list in s2_components:
        component_prefix = component_info["prefix"]
        component_name = component_info["name"]
        
        # This part of the code generates the terms for the final equation printout
        expanded_terms_str = []
        for term in term_list:
            final_term_literals = component_prefix + term
            num_literals = len(final_term_literals)
            multiplications = num_literals - 1
            total_multiplications += multiplications
            multiplication_counts.append(multiplications)
            expanded_terms_str.append(" * ".join(final_term_literals))
        
        final_s2_expression_parts.append(" + ".join(expanded_terms_str))

    print("s2 = " + " + \n     ".join(final_s2_expression_parts))
    print("\nCounting the multiplications for each product term (a term with k literals has k-1 multiplications):")
    
    # Print the sum of all multiplication counts
    sum_str = " + ".join(map(str, multiplication_counts))
    print(f"Total Multiplications = {sum_str} = {total_multiplications}")

solve()