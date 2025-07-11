from sympy import symbols
from sympy.logic.boolalg import And, Or, Xor, to_dnf, Not

def solve_s2_multiplications():
    """
    This function calculates the number of multiplication operations in the
    fully expanded expression for the sum bit s2 of a 3-bit binary adder.
    """
    # Define the boolean variables for the two 3-bit numbers A and B
    a0, a1, a2, b0, b1, b2 = symbols('a0, a1, a2, b0, b1, b2')

    # c1 is the carry-out from the LSB (position 0) addition
    c1 = And(a0, b0)

    # c2 is the carry-out from the position 1 addition.
    # It can be expressed as a majority function of its inputs: a1, b1, c1
    c2 = Or(And(a1, b1), And(a1, c1), And(b1, c1))

    # s2 is the sum bit for position 2, which is the XOR of its inputs a2, b2, and c2
    # In SOP form: s2 = (a2'b2'c2) + (a2'b2c2') + (a2b2'c2') + (a2b2c2)
    s2 = Xor(a2, b2, c2)

    # Convert the final expression for s2 into its Disjunctive Normal Form (DNF),
    # which is another name for a fully expanded Sum-of-Products (SOP) expression.
    s2_dnf = to_dnf(s2, simplify=True)
    
    # Helper function to format a literal for printing (e.g., Not(a0) -> a0')
    def format_literal(lit):
        s = str(lit)
        if s.startswith('~'):
            return f"{s[1:]}'"
        return s

    total_multiplications = 0
    
    # The result of to_dnf is an Or of Ands. We extract the And terms.
    if isinstance(s2_dnf, Or):
        product_terms = s2_dnf.args
    else: # Handle cases where the DNF is a single term
        product_terms = [s2_dnf]
    
    equation_parts = []

    # Iterate through each product term in the SOP expression
    for term in product_terms:
        # Each term is usually an And object. If not, it's a single literal.
        if isinstance(term, And):
            literals = term.args
            num_literals = len(literals)
            # A product of n literals requires n-1 multiplication operations
            if num_literals > 1:
                total_multiplications += (num_literals - 1)
            
            # Format the term for display, e.g., "a0' * a1 * b2"
            formatted_literals = sorted([format_literal(l) for l in literals])
            equation_parts.append(" * ".join(formatted_literals))

        else: # A single literal term has 0 multiplications.
            equation_parts.append(format_literal(term))
    
    # Sort the final terms alphabetically for a clean final equation
    equation_parts.sort()
    final_equation = " + ".join(equation_parts)

    print("The fully expanded expression for s2 is:")
    print(f"s2 = {final_equation}")
    print("\nThis equation describes the logic for the s2 bit. To find the total number of")
    print("multiplication operations, we count the number of literals in each product term")
    print("and subtract one. For example, a term 'a*b*c' has 2 multiplications.")
    print("\nTotal number of multiplication operations in the final equation:")
    print(total_multiplications)

solve_s2_multiplications()