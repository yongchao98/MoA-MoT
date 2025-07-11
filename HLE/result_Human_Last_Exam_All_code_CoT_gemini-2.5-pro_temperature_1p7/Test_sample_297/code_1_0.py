import sympy
from sympy import symbols
from sympy.logic.boolalg import Or, And, Not

def solve_s2_expression():
    """
    This function derives the fully expanded expression for the sum bit s2 of a
    3-bit binary adder and counts the number of multiplication operations.
    """
    # Define the binary digits of A and B as symbolic variables
    # a0, b0 are the least significant bits.
    a0, a1, a2 = symbols('a0, a1, a2')
    b0, b1, b2 = symbols('b0, b1, b2')

    # Step 1: Define the carry-out from the LSB stage (c1)
    # c1 = a0 AND b0
    c1 = And(a0, b0)

    # Step 2: Define the carry-out from the middle stage (c2)
    # c2 = (a1 AND b1) OR (a1 AND c1) OR (b1 AND c1)
    # sympy automatically substitutes c1, resulting in:
    # c2 = (a1 & b1) | (a0 & b0 & a1) | (a0 & b0 & b1)
    # This is the minimal Sum-of-Products (SOP) form for c2.
    c2_sop = Or(And(a1, b1), And(a1, c1), And(b1, c1))
    
    # Step 3: Define the complement of c2 (c2'), also in minimal SOP form.
    # This is needed for the XOR expansion.
    c2_prime_sop = sympy.logic.to_dnf(Not(c2_sop), simplify=True)

    # Step 4: Define s2 using the SOP expansion of a 3-input XOR gate:
    # s2 = a2 XOR b2 XOR c2 expands to:
    # (a2'b2'c2) + (a2'b2c2') + (a2b2'c2') + (a2b2c2)
    # We construct this by ANDing the terms with the SOP forms of c2 and c2'.
    term1 = And(Not(a2), Not(b2), c2_sop)
    term2 = And(Not(a2), b2, c2_prime_sop)
    term3 = And(a2, Not(b2), c2_prime_sop)
    term4 = And(a2, b2, c2_sop)
    
    # Combine the four major components into a single expression for s2
    s2_expression = Or(term1, term2, term3, term4)

    # Step 5: "Fully expand" the expression by distributing all terms using sympy.expand().
    # This results in the final Sum-of-Products (SOP) form with no nested logic.
    s2_expanded_sop = sympy.expand(s2_expression)
    
    # Get the list of individual product terms from the final SOP expression.
    if isinstance(s2_expanded_sop, sympy.logic.boolalg.Or):
        product_terms = s2_expanded_sop.args
    else: # In case the expression simplifies to a single term
        product_terms = [s2_expanded_sop]
    
    total_multiplications = 0
    
    print("The final, fully expanded expression for s2 is the sum (OR) of the following product (AND) terms:")
    print("Note: In boolean algebra, 'multiplication' is the AND operation, and complement is shown with '~'.")
    print("-" * 80)
    
    # Step 6: Loop through each product term to count the multiplications.
    # We sort the terms to ensure the output is in a consistent, readable order.
    for term in sorted(product_terms, key=str):
        # A product term can be a single literal or an And(...) object in sympy.
        if isinstance(term, sympy.logic.boolalg.And):
            num_literals = len(term.args)
        else:
            num_literals = 1
        
        # A product of 'k' literals (e.g., a&b&c) requires 'k-1' multiplication operations.
        multiplications_in_term = max(0, num_literals - 1)
        total_multiplications += multiplications_in_term
        
        # Print the breakdown for each term.
        print(f"Term: {str(term):<35} => Literals: {num_literals}, Multiplications: {multiplications_in_term}")

    print("-" * 80)
    print(f"The total number of multiplication operations is: {total_multiplications}")

if __name__ == '__main__':
    solve_s2_expression()