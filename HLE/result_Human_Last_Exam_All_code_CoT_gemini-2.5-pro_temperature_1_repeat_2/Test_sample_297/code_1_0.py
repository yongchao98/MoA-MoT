import sympy
from sympy.logic.boolalg import And, Or, Not, to_sop

def solve_and_count_multiplications():
    """
    This function programmatically derives the fully expanded expression for the 
    sum bit s2 of a 3-bit binary adder and counts the total number of 
    multiplication (AND) operations.
    """

    # Plan Step 1: Model the 3-bit Adder by defining the input bits.
    # A = a2a1a0, B = b2b1b0
    a0, a1, a2, b0, b1, b2 = sympy.symbols('a0, a1, a2, b0, b1, b2')

    # Plan Step 2: Derive the expression for s2 step-by-step.
    
    # The carry-out from the first stage (LSB: a0 + b0) is c1.
    c1 = a0 & b0

    # The carry-out from the second stage (a1 + b1 + c1) is c2.
    # The formula for carry-out of a full adder is (x&y) | (x&z) | (y&z).
    c2_expression = (a1 & b1) | (a1 & c1) | (b1 & c1)

    # The sum bit s2 is the sum output from the third stage (a2 + b2 + c2).
    # The formula for the sum of a full adder is x ^ y ^ z (XOR).
    # The SOP form of x^y^z is (x'y'z) + (x'yz') + (xy'z') + (xyz).
    # We build this expression for s2 using a2, b2, and the derived c2.
    s2_full_expression = (
        (~a2 & ~b2 & c2_expression) |
        (~a2 & b2 & ~c2_expression) |
        (a2 & ~b2 & ~c2_expression) |
        (a2 & b2 & c2_expression)
    )

    # Plan Step 3: Convert the logical expression for s2 into its fully expanded
    # Sum of Products (SOP) form. This is the "fully expanded expression"
    # mentioned in the problem.
    s2_sop = to_sop(s2_full_expression, force=True)

    # Plan Step 4: Count the multiplications in the final expression.
    # Each AND operation is a multiplication. A product term with N variables
    # (e.g., a*b*c) contains N-1 multiplication operations.
    
    total_multiplications = 0
    
    # The SOP form is an Or(...) of And(...) terms. We extract the terms.
    if isinstance(s2_sop, Or):
        product_terms = s2_sop.args
    elif isinstance(s2_sop, And):
        product_terms = [s2_sop]
    else:
        product_terms = []

    print("The final equation for s2 is a sum of the following product terms:")
    for term in product_terms:
        # sympy uses `&` for AND, `|` for OR, `~` for NOT.
        # This corresponds to multiplication, addition, and complement.
        print(f"  {term}")
    
    # Create a list to hold the multiplication count for each term.
    term_mult_counts = []
    for term in product_terms:
        if isinstance(term, And):
            num_literals = len(term.args)
        else: # A term with a single variable like `a0`
            num_literals = 1
        
        # A term with N literals has N-1 multiplications.
        mult_count = num_literals - 1
        if mult_count > 0:
            total_multiplications += mult_count
            term_mult_counts.append(str(mult_count))

    calculation_str = " + ".join(term_mult_counts)
    
    print("\nEach number in the final equation below represents the count of multiplications for a single term in the expanded expression for s2.")
    print(f"Final Equation for Total Multiplications: {calculation_str}")
    print(f"\nThe total number of multiplication operations is: {total_multiplications}")

# Run the solver function
solve_and_count_multiplications()