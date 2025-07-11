import sympy
from sympy.logic.boolalg import Or, And, Not
from sympy.abc import A, B, C

def main():
    """
    Calculates the number of multiplications in the fully expanded
    boolean expression for the sum bit s2 of a 3-bit adder.
    """

    # Define symbols for the binary digits. Using a trailing underscore for complement.
    # e.g., a2_ represents the complement of a2.
    a0, a1, a2, b0, b1, b2 = sympy.symbols('a0, a1, a2, b0, b1, b2')
    a0_, a1_, a2_, b0_, b1_, b2_ = sympy.symbols("a0', a1', a2', b0', b1', b2'")

    # Mapping from complement symbol to original for string replacement
    complement_map = {
        a0_: "a0'", a1_: "a1'", a2_: "a2'",
        b0_: "b0'", b1_: "b1'", b2_: "b2'"
    }

    # Define the carry expressions based on boolean logic for an adder
    # c1 is the carry-out from the 0th bit (a0 + b0)
    c1 = a0 * b0

    # c2 is the carry-out from the 1st bit (a1 + b1 + c1)
    c2 = sympy.expand((a1 * b1) + (a1 * c1) + (b1 * c1))

    # To get the complement c2', we first need c1'
    c1_ = a0_ + b0_ # De Morgan's Law on c1 = a0*b0
    
    # The complement of a carry-out is c_out' = a'*b' + a'*c_in' + b'*c_in'
    c2_ = sympy.expand((a1_ * b1_) + (a1_ * c1_) + (b1_ * c1_))
    
    # s2 is the sum bit from the 2nd stage (a2 + b2 + c2)
    # The sum bit is the XOR of the inputs: s2 = a2 XOR b2 XOR c2
    # This expands to: s2 = a2'b2'c2 + a2'b2c2' + a2b2'c2' + a2b2c2
    s2_expr = (a2_ * b2_ * c2) + (a2_ * b2 * c2_) + \
              (a2 * b2_ * c2_) + (a2 * b2 * c2)

    # Fully expand the expression for s2 into its sum-of-products form
    s2_expanded = sympy.expand(s2_expr)

    # Format the final expression for printing
    # Replace complement symbols (e.g., a2_) with a prime (e.g., a2')
    s2_str = str(s2_expanded)
    for comp_var, prime_var in complement_map.items():
        s2_str = s2_str.replace(str(comp_var), prime_var)
    s2_str = s2_str.replace("*", " * ")

    print("The fully expanded expression for s2 is:")
    print("s2 = " + s2_str)
    print("\n" + "="*50)
    print("Calculation of multiplication operations for each term:")
    print("="*50)

    total_multiplications = 0
    mult_counts = []
    
    # Iterate through each product term in the expanded sum
    # .as_ordered_terms() splits the sum into its component terms
    for term in s2_expanded.as_ordered_terms():
        # .as_ordered_factors() splits a product into its variables
        factors = term.as_ordered_factors()
        num_vars = len(factors)
        
        # A product of n variables requires n-1 multiplications
        num_muls = 0
        if num_vars > 1:
            num_muls = num_vars - 1
        
        total_multiplications += num_muls
        mult_counts.append(str(num_muls))

        term_str = str(term)
        for comp_var, prime_var in complement_map.items():
            term_str = term_str.replace(str(comp_var), prime_var)
        term_str = term_str.replace("*", "")
        print(f"Term '{term_str}': {num_vars} variables => {num_muls} multiplications")

    sum_str = " + ".join(mult_counts)
    print("\n" + "="*50)
    print(f"Total multiplications = {sum_str} = {total_multiplications}")
    print("="*50)

if __name__ == "__main__":
    main()
<<<52>>>