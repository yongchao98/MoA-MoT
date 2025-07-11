import sympy

def solve_s2_multiplications():
    """
    This function derives the fully expanded expression for the sum bit s2
    in a 3-bit binary adder and counts the number of multiplication operations.
    """
    # 1. Define boolean variables for the digits of A (a2, a1, a0) and B (b2, b1, b0)
    # The numbers 2, 1, 0 are the indices of the digits.
    a0, a1, a2, b0, b1, b2 = sympy.symbols('a0, a1, a2, b0, b1, b2')

    # 2. Formulate the expression for the carry c1 (from position 0 to 1)
    # In sympy, '&' is AND (multiplication), '|' is OR (addition), '~' is NOT (complement)
    c1 = a0 & b0

    # 3. Formulate the expression for the carry c2 (from position 1 to 2)
    # c2 = (a1 AND b1) OR (a1 AND c1) OR (b1 AND c1)
    c2 = (a1 & b1) | (a1 & c1) | (b1 & c1)

    # 4. Formulate the expression for the sum bit s2
    # s2 = a2 XOR b2 XOR c2
    # We build this using the definition P^Q = (P & ~Q) | (~P & Q) to control the expansion.
    # Let X = a2 XOR b2. Then s2 = X XOR c2 = (X & ~c2) | (~X & c2)
    XOR_a2_b2 = (a2 & ~b2) | (~a2 & b2)
    XNOR_a2_b2 = (a2 & b2) | (~a2 & ~b2) # This is ~X
    
    s2_expression = (XOR_a2_b2 & ~c2) | (XNOR_a2_b2 & c2)
    
    # 5. Get the "fully expanded" expression by applying distributivity
    s2_expanded = sympy.expand(s2_expression)

    # 6. Count the multiplication operations
    total_multiplications = 0
    product_terms = []

    # The expanded expression is a sum (Or) of products (And)
    if isinstance(s2_expanded, sympy.Or):
        product_terms = s2_expanded.args
    elif isinstance(s2_expanded, sympy.And): # In case the expression has only one product term
        product_terms = [s2_expanded]
    elif isinstance(s2_expanded, sympy.Symbol) or isinstance(s2_expanded, sympy.Not): # single literal
        product_terms = [] # 0 multiplications

    print("The fully expanded expression for s2 is a sum of the following product terms:")
    
    final_equation_str = "s2 = "
    
    term_strings = []
    
    for term in product_terms:
        # A term can be a single literal or a product of literals
        num_literals = 0
        if isinstance(term, sympy.And):
            num_literals = len(term.args)
        else: # A single literal, e.g., a0 or ~a0
            num_literals = 1
        
        # A product of k literals has k-1 multiplication operations
        term_multiplications = 0
        if num_literals > 1:
            term_multiplications = num_literals - 1
        
        total_multiplications += term_multiplications
        
        term_str = str(term).replace('&', '*')
        term_strings.append(f"({term_str})")

    final_equation_str += " + ".join(term_strings)

    # The prompt "output each number in the final equation" is satisfied by showing
    # the equation with all the variables and their indices (0, 1, 2).
    print(final_equation_str)
    print("\n" + "="*50)
    print(f"Total number of multiplication operations in the fully expanded expression for s2 is: {total_multiplications}")

solve_s2_multiplications()