import sympy
from sympy.logic.boolalg import to_dnf, And, Or, Not

def solve_and_print_s2_expression():
    """
    This function derives the expression for the sum bit s2 of a 3-bit binary adder,
    expands it to its full Sum-of-Products form, and counts the multiplication operations.
    """
    # Define the binary digits of A and B as symbolic variables
    # a2, a1, a0 are the digits for number A
    # b2, b1, b0 are the digits for number B
    a2, a1, a0 = sympy.symbols('a2 a1 a0')
    b2, b1, b0 = sympy.symbols('b2 b1 b0')

    # Step 1: Define the carry-out c1 from the first bit (a0, b0)
    # c1 = a0 AND b0
    c1 = And(a0, b0)

    # Step 2: Define the carry-out c2 from the second bit (a1, b1)
    # c2 depends on a1, b1, and the carry-in c1
    # c2 = (a1 AND b1) OR (a1 AND c1) OR (b1 AND c1)
    c2_expr_with_c1 = Or(And(a1, b1), And(a1, c1), And(b1, c1))

    # Substitute c1 into the expression for c2 to make it explicit
    # in terms of the base digits.
    c2 = c2_expr_with_c1.subs({c1: And(a0, b0)})

    # Step 3: Define the sum bit s2
    # s2 is the 3-input XOR of a2, b2, and the carry-in c2
    s2_expr = a2 ^ b2 ^ c2

    # Step 4: Convert the expression for s2 into Sum of Products form (DNF)
    # The to_dnf function with simplify=True gives the minimal SOP form.
    s2_dnf = to_dnf(s2_expr, simplify=True)

    # Step 5: Count the total multiplications and build the equation string
    # A product term with N literals requires N-1 multiplications.
    total_multiplications = 0
    equation_parts = []
    
    product_terms = []
    if isinstance(s2_dnf, Or):
        product_terms = s2_dnf.args
    else: # Handles case where DNF is a single term
        product_terms = [s2_dnf]
    
    # Process each product term
    for term in product_terms:
        literals = []
        if isinstance(term, And):
            literals = term.args
        else:
            literals = [term]
        
        # Build the string for the current term for the final equation
        term_str_parts = []
        for literal in sorted(literals, key=str): # Sort for consistent output
            if isinstance(literal, Not):
                # Represent complement with a prime (') as in the problem description
                term_str_parts.append(f"{literal.args[0]}'")
            else:
                term_str_parts.append(str(literal))
        equation_parts.append(" * ".join(term_str_parts))

    # Sort the final equation terms for readability
    equation_parts.sort()
    final_equation = "s2 = " + " + ".join(equation_parts)

    print("The fully expanded expression for s2 is:")
    print(final_equation)
    print("\n" + "="*80)
    print("Counting the multiplication operations for each product term:")
    print("A term with N literals (variables) requires N-1 multiplications.")
    print("="*80)

    # Now, iterate through the generated string terms to count and display
    for term_str in equation_parts:
        # The number of literals is one more than the number of '*' separators
        num_literals = term_str.count('*') + 1
        if num_literals > 1:
            term_mults = num_literals - 1
            total_multiplications += term_mults
            print(f"- Term '{term_str}': {num_literals} literals => {term_mults} multiplications.")
        else:
            # A term with 1 literal has 0 multiplications
            print(f"- Term '{term_str}': {num_literals} literal => 0 multiplications.")

    print("="*80)
    print(f"\nThe total number of multiplication operations is the sum of the multiplications from each term.")
    print(f"Total Multiplications = {total_multiplications}")


if __name__ == '__main__':
    solve_and_print_s2_expression()