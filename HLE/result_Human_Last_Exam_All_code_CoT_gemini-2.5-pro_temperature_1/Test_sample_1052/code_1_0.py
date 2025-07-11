import itertools

def zhigalkin_polynomial(a, b, c, d):
    """
    Calculates the value of the given Zhigalkin polynomial.
    P(a,b,c,d) = c ⊕ d ⊕ (b∧c) ⊕ (a∧d) ⊕ (a∧c) ⊕ (a∧b∧d) ⊕ (a∧b∧c)
    In Python, bitwise XOR is ^ and bitwise AND is &.
    """
    return (c ^ d ^ (b & c) ^ (a & d) ^ (a & c) ^ (a & b & d) ^ (a & b & c))

def derived_formula(a, b, c, d):
    """
    Calculates the value of the derived Boolean formula.
    Formula: (b ∧ d) ∨ ((a ↓ b) ∧ ¬(c ↔ d))
    This is implemented using Python's bitwise operators for 0/1 logic.
    """
    # Define logical operators using 0/1 logic
    op_not = lambda x: 1 - x
    op_or = lambda x, y: x | y
    op_nor = lambda x, y: 1 - (x | y)  # a ↓ b is equivalent to ¬(a ∨ b)
    op_iff = lambda x, y: 1 - (x ^ y)  # a ↔ d is ¬(a ⊕ b)
    op_and = lambda x, y: x & y
    
    # Break down the formula: (b ∧ d) ∨ ((a ↓ b) ∧ ¬(c ↔ d))
    term1 = op_and(b, d)
    
    term2_part1 = op_nor(a, b)
    term2_part2 = op_not(op_iff(c, d))
    term2 = op_and(term2_part1, term2_part2)
    
    result = op_or(term1, term2)
    return result

def main():
    """
    Main function to verify the equivalence and print the result.
    """
    print("Verifying the equivalence of the polynomial and the derived formula.")
    print("-" * 50)
    print("a b c d | Polynomial | Formula | Match?")
    print("-" * 50)
    
    all_match = True
    # Iterate through all 16 possible boolean inputs
    for a, b, c, d in itertools.product([0, 1], repeat=4):
        poly_result = zhigalkin_polynomial(a, b, c, d)
        formula_result = derived_formula(a, b, c, d)
        match = "Yes" if poly_result == formula_result else "No"
        if poly_result != formula_result:
            all_match = False
        print(f"{a} {b} {c} {d} |      {poly_result}     |    {formula_result}    |  {match}")

    print("-" * 50)
    if all_match:
        print("Verification successful: The polynomial and the formula are equivalent.")
    else:
        print("Verification failed: There is a mismatch.")

    # To express the formula (b ∧ d) ∨ ((a ↓ b) ∧ ¬(c ↔ d)) using the allowed operators
    # (¬, ↑, ↓, ↔, →, ∨), we replace ∧ with ¬(A ↑ B).
    # (b ∧ d) becomes ¬(b ↑ d)
    # ((a ↓ b) ∧ ¬(c ↔ d)) becomes ¬((a ↓ b) ↑ ¬(c ↔ d))
    final_formula_str = "¬(b ↑ d) ∨ ¬((a ↓ b) ↑ ¬(c ↔︎ d))"
    
    print("\nThe Boolean formula derived from the polynomial is:")
    print(final_formula_str)

if __name__ == "__main__":
    main()
