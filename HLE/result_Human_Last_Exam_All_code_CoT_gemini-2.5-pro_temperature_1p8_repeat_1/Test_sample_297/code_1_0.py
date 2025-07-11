import sympy
from sympy.logic import SOPform
from sympy.logic.boolalg import And, Or, Xor, Not

def solve_s2_multiplications():
    """
    This function calculates the number of multiplications in the fully expanded
    Sum-of-Products (SOP) expression for the sum bit s2 of a 3-bit binary adder.
    """
    # Define the boolean variables for the two 3-bit numbers A and B
    a0, a1, a2, b0, b1, b2 = sympy.symbols('a0, a1, a2, b0, b1, b2')

    # Step 1: Define the carry-in to bit 2, which is the carry-out from bit 1.
    # c1 = a0 AND b0
    c1 = And(a0, b0)
    # c2 = (a1 AND b1) OR (c1 AND (a1 XOR b1))
    c2_expr = Or(And(a1, b1), And(c1, Xor(a1, b1)))

    # Step 2: Define the sum bit s2.
    # s2 = a2 XOR b2 XOR c2
    s2_expr = Xor(a2, Xor(b2, c2_expr))

    # Step 3: Convert the expression for s2 into its minimal Sum-of-Products form.
    # This automatically expands all logic and simplifies to the final expression.
    s2_sop = SOPform.from_formula(s2_expr, simplify=True)
    
    # Step 4: Iterate through the final product terms and count the multiplications.
    # A product term with 'n' literals requires 'n-1' multiplications.
    total_multiplications = 0
    product_terms = []
    
    if isinstance(s2_sop, Or):
        product_terms = s2_sop.args
    elif isinstance(s2_sop, And):
        product_terms = [s2_sop]

    # This dictionary will hold the count of terms for each literal size.
    # e.g., {4: 6, 5: 8} means 6 terms with 4 literals, 8 terms with 5 literals.
    term_counts_by_literals = {}

    for term in product_terms:
        num_literals = len(term.args)
        term_counts_by_literals[num_literals] = term_counts_by_literals.get(num_literals, 0) + 1

    print("Counting multiplications in the fully expanded expression for s2:")
    
    calc_str_parts = []
    # Sort by number of literals for a clear breakdown
    for num_lit, count in sorted(term_counts_by_literals.items()):
        mults_for_group = count * (num_lit - 1)
        total_multiplications += mults_for_group
        
        # This part fulfills the requirement: "output each number in the final equation"
        print(f"There are {count} terms with {num_lit} literals each.")
        print(f"Contribution: {count} * ({num_lit} - 1) = {mults_for_group} multiplications.")
        calc_str_parts.append(str(mults_for_group))
    
    print("\nFinal Calculation: " + " + ".join(calc_str_parts) + f" = {total_multiplications}")
    
    print("\nTherefore, the total number of multiplication operations is:")
    print(total_multiplications)

if __name__ == '__main__':
    solve_s2_multiplications()