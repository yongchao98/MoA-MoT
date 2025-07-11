import sympy
from sympy import symbols
from sympy.logic import simplify_logic
from sympy.printing.str import StrPrinter

def main():
    """
    Calculates the number of multiplications in the fully expanded 
    boolean expression for the s2 bit of a 3-bit adder.
    """
    # 1. Define the boolean variables for the two 3-bit numbers A and B.
    # a0, a1, a2 are bits for A; b0, b1, b2 are bits for B.
    a0, a1, a2, b0, b1, b2 = symbols('a0,a1,a2, b0,b1,b2')
    
    # 2. Define the carry logic based on a standard full adder.
    # c1 is the carry-out from adding a0 and b0.
    c1 = a0 & b0
    
    # c2 is the carry-out from adding a1, b1, and the carry-in c1.
    c2 = (a1 & b1) | (a1 & c1) | (b1 & c1)
    
    # 3. Define the sum bit s2.
    # s2 is the result of XORing a2, b2, and the carry-in c2.
    s2_expr = a2 ^ b2 ^ c2
    
    # 4. Simplify the expression for s2 into its minimal Sum-of-Products form (DNF).
    # This is the "fully expanded expression" using only OR (addition) and AND (multiplication).
    s2_dnf = simplify_logic(s2_expr, form='dnf')

    # To format with overlines for complements, we'll create a custom printer
    class OverlinePrinter(StrPrinter):
        def _print_Not(self, expr):
            # args[0] is the argument of Not
            return f"{self._print(expr.args[0])}'"

    custom_printer = OverlinePrinter()
    
    print("The fully expanded expression for s2 is a sum of the following product terms:")
    
    total_mults = 0
    # The result s2_dnf is an Or(...) of several And(...) terms.
    # We retrieve the list of And terms from the arguments of the Or.
    product_terms = s2_dnf.args if isinstance(s2_dnf, sympy.logic.boolalg.Or) else [s2_dnf]

    # 5. Iterate through each product term to count literals and multiplications.
    equation_parts = []
    for i, term in enumerate(product_terms):
        num_literals = len(term.args) if isinstance(term, sympy.logic.boolalg.And) else 1
        
        # A product of k literals has k-1 multiplication operations.
        mults_in_term = 0
        if num_literals > 1:
            mults_in_term = num_literals - 1
            
        total_mults += mults_in_term
        
        # Format the term with complements shown as an overline (e.g., a0')
        term_str = "".join(custom_printer.doprint(arg) for arg in (term.args if isinstance(term, sympy.logic.boolalg.And) else [term]))
        equation_parts.append(term_str)
        
        # Per the prompt, "output each number in the final equation"
        print(f"- Term {i+1}: {term_str} (has {num_literals} literals, requires {mults_in_term} multiplications)")
        
    print("\nThe final equation is s2 = " + " + ".join(equation_parts))
    print(f"\nThe total number of multiplication operations is the sum of the multiplications from each term.")
    print(f"Total multiplications = {total_mults}")

if __name__ == '__main__':
    main()
<<<38>>>