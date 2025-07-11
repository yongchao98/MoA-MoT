import sympy
from sympy.logic.boolalg import And, Or, Not, Xor, to_sop
from sympy import symbols

def solve_s2_multiplications():
    """
    This function calculates the number of multiplications in the fully expanded
    Sum-of-Products (SOP) expression for the sum bit s2 of a 3-bit binary adder.
    """
    # Step 1: Define the binary digits as symbols.
    # A = a2,a1,a0 and B = b2,b1,b0
    a0, a1, a2 = symbols('a0, a1, a2')
    b0, b1, b2 = symbols('b0, b1, b2')

    # Step 2: Define the logic for the carry bits.
    # The initial carry-in c0 is assumed to be 0.
    c1 = And(a0, b0)
    # The carry-in to the third stage (c2) depends on a1, b1, and c1.
    c2 = Or(And(a1, b1), And(a1, c1), And(b1, c1))
    
    # Fully expand c2 in terms of the primary inputs.
    # sympy.logic.simplify_logic can be used, but subs is more direct here.
    c2_expanded = c2.subs(symbols('c1'), c1)

    # Step 3: Define the sum bit s2.
    # s2 is the XOR of a2, b2, and the carry-in c2.
    s2 = Xor(a2, b2, c2_expanded)

    # Step 4: Convert the expression for s2 to its fully expanded SOP form.
    # The `force=True` argument is needed to expand XOR/XNOR expressions.
    s2_sop = to_sop(s2, force=True)

    # Step 5: Count the multiplication operations.
    # The result s2_sop is an Or object, where each argument is a product term.
    # We check if the result is an Or, otherwise it's a single term.
    if isinstance(s2_sop, Or):
        product_terms = s2_sop.args
    else:
        product_terms = [s2_sop]

    mult_counts = []
    total_multiplications = 0

    for term in product_terms:
        # A product term is an And object. The number of literals is its number of arguments.
        # If the term is not an And (i.e., a single literal), it has 1 literal.
        if isinstance(term, And):
            num_literals = len(term.args)
        else:
            num_literals = 1
        
        # A product of 'n' literals requires 'n-1' multiplication operations.
        # A term with a single literal requires 0 multiplications.
        if num_literals > 1:
            multiplications_in_term = num_literals - 1
            mult_counts.append(multiplications_in_term)
            total_multiplications += multiplications_in_term

    # Sort the list of counts for a consistent and readable output.
    mult_counts.sort(reverse=True)

    # Step 6: Print the results as per the plan.
    print("The final expression for s2 is a sum of multiple product terms.")
    print("A product term with 'n' variables (e.g., a2*b1'*a0) requires 'n-1' multiplications.")
    print("\nThe number of multiplications for each of the 16 product terms in the expanded expression are:")
    print(mult_counts)
    
    # The prompt asks to output each number in the final equation.
    equation_str = " + ".join(map(str, mult_counts))
    print("\nThe total number of multiplications is the sum of these individual counts:")
    print(f"Total Multiplications = {equation_str} = {total_multiplications}")

    return total_multiplications

# Execute the function and print the final answer.
final_answer = solve_s2_multiplications()
print(f"\n<<<{final_answer}>>>")
