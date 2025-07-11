from sympy import symbols
from sympy.logic.boolalg import to_sop, And, Or, Not, Xor

def solve_s2_expression():
    """
    This function derives the fully expanded expression for the s2 bit
    of a 3-bit binary adder and counts the multiplication operations.
    """
    # Define the input binary digits as symbols
    a0, a1, a2, b0, b1, b2 = symbols('a0, a1, a2, b0, b1, b2')

    # Define the logic for the 3-bit adder
    # c1 is the carry-out from the first bit (LSB) addition
    c1 = And(a0, b0)

    # c2 is the carry-out from the second bit addition
    # c2 = (a1 AND b1) OR (a1 AND c1) OR (b1 AND c1)
    c2 = Or(And(a1, b1), And(a1, c1), And(b1, c1))

    # s2 is the sum bit from the third bit addition
    # s2 = a2 XOR b2 XOR c2
    s2_expr = Xor(a2, b2, c2)

    # Convert the expression for s2 to its minimal Sum of Products (SOP) form.
    # This represents the "fully expanded expression".
    s2_sop = to_sop(s2_expr, force=True)

    total_multiplications = 0
    
    # The result of to_sop can be a single term or an Or clause.
    product_terms = s2_sop.args if isinstance(s2_sop, Or) else [s2_sop]
    
    # We will build the string for the final equation and count multiplications
    expression_parts = []
    calculation_details = []

    for term in product_terms:
        # Determine the number of literals in the product term
        if isinstance(term, And):
            num_literals = len(term.args)
            literals_str = [str(lit).replace("~", "") + "'" if isinstance(lit, Not) else str(lit) for lit in term.args]
        elif isinstance(term, Not) or isinstance(term, symbols.Symbol):
            num_literals = 1
            lit_str = str(term)
            literals_str = [lit_str.replace("~", "") + "'" if isinstance(term, Not) else lit_str]
        else:
            # This case should not be reached for a valid SOP form
            continue
        
        # Add to the full equation string
        expression_parts.append(" ⋅ ".join(literals_str))
        
        # A product of k literals requires k-1 multiplications
        muls = num_literals - 1
        if muls > 0:
            total_multiplications += muls
        
        calculation_details.append(f"  - Term '{' ⋅ '.join(literals_str)}': {num_literals} literals -> {muls} multiplications")

    print("The fully expanded expression for s2 (using ⋅ for multiplication, + for addition, and ' for complement) is:")
    print("s2 = " + " + ".join(expression_parts))
    print("\nEach term in the sum is a product. The number of multiplications for a product of k digits is k-1.")
    print("The terms in the expanded expression and their corresponding multiplication counts are:")
    print("\n".join(calculation_details))
    print(f"\nTotal number of multiplication operations = {total_multiplications}")

solve_s2_expression()