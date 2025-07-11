def solve_s2_multiplications():
    """
    Calculates the number of multiplication operations in the fully expanded
    boolean expression for the sum bit s2 of a 3-bit binary adder.
    """

    # Step 4 & 5: Define the terms in the minimal SOP for c2 and c2'
    # c2 = a1b1 + a1a0b0 + b1a0b0
    c2_terms = [
        {'a1', 'b1'},
        {'a1', 'a0', 'b0'},
        {'b1', 'a0', 'b0'}
    ]

    # c2' = a1'b1' + a1'a0' + a1'b0' + b1'a0' + b1'b0'
    c2_prime_terms = [
        {"a1'", "b1'"},
        {"a1'", "a0'"},
        {"a1'", "b0'"},
        {"b1'", "a0'"},
        {"b1'", "b0'"}
    ]

    # Step 1: The s2 expression is composed of four main parts
    # s2 = (a2'b2')c2 + (a2'b2)c2' + (a2b2')c2' + (a2b2)c2
    s2_parts = [
        ({"a2'", "b2'"}, c2_terms),      # (a2'b2')c2
        ({"a2'", "b2"}, c2_prime_terms), # (a2'b2)c2'
        ({"a2", "b2'"}, c2_prime_terms), # (a2b2')c2'
        ({"a2", "b2"}, c2_terms)       # (a2b2)c2
    ]

    total_multiplications = 0
    final_expression_terms = []

    # Step 6 & 7: Expand s2 and count multiplications
    print("Deriving the fully expanded expression for s2 and counting multiplications:\n")

    for prefix_lits, term_list in s2_parts:
        # Example: prefix_lits = {"a2'", "b2'"}, term_list = c2_terms
        for term_lits in term_list:
            # Combine the prefix with each term
            final_term = prefix_lits.union(term_lits)
            final_expression_terms.append(final_term)
            
            # A product term with N literals has N-1 multiplications
            num_literals = len(final_term)
            multiplications = num_literals - 1
            total_multiplications += multiplications
            
            # Print the breakdown for clarity
            term_str = " * ".join(sorted(list(final_term)))
            print(f"Term: {term_str}")
            print(f"  Literals: {num_literals}, Multiplications: {multiplications}\n")

    print("-----------------------------------------")
    print(f"Total number of product terms in the final expression for s2: {len(final_expression_terms)}")
    print(f"Total number of multiplication operations: {total_multiplications}")

solve_s2_multiplications()
<<<52>>>