def zhigalkin_polynomial(a, b, c, d):
    """
    Calculates the value of the given Zhigalkin polynomial.
    In Python, a^b is XOR, a&b is AND.
    The polynomial is:
    d^c^(b&c)^(a&d)^(a&c)^(a&b&d)^(a&b&c)
    """
    # Convert boolean True/False to 1/0 if necessary
    a, b, c, d = int(a), int(b), int(c), int(d)
    
    term1 = d
    term2 = c
    term3 = b & c
    term4 = a & d
    term5 = a & c
    term6 = a & b & d
    term7 = a & b & c
    
    return term1 ^ term2 ^ term3 ^ term4 ^ term5 ^ term6 ^ term7

def derived_boolean_formula(a, b, c, d):
    """
    Calculates the value of the derived Boolean formula:
    ((a ∨ b) ∨ (c ↔ d)) → ¬(b → ¬d)
    """
    # Python's implementation of logical operators:
    # a or b  (OR, ∨)
    # a == b  (IFF, ↔)
    # not a   (NOT, ¬)
    # (not a) or b (IMPLIES, →)
    
    # Evaluate sub-expressions
    a_or_b = a or b
    c_iff_d = (c == d)
    
    # Antecedent of the main implication
    antecedent = a_or_b or c_iff_d
    
    # Consequent of the main implication: ¬(b → ¬d)
    # b → ¬d is equivalent to (¬b ∨ ¬d)
    b_implies_not_d = (not b) or (not d)
    consequent = not b_implies_not_d
    
    # Final implication: antecedent → consequent
    # This is equivalent to ¬antecedent ∨ consequent
    return (not antecedent) or consequent

def solve():
    """
    Verifies the equivalence of the polynomial and the formula and prints the result.
    """
    print("Verifying the equivalence of the polynomial and the derived formula...")
    print("-" * 60)
    print("a b c d | Poly Output | Formula Output | Match?")
    print("-" * 60)
    
    all_match = True
    for a in [False, True]:
        for b in [False, True]:
            for c in [False, True]:
                for d in [False, True]:
                    poly_result = zhigalkin_polynomial(a, b, c, d)
                    formula_result = derived_boolean_formula(a, b, c, d)
                    
                    match_str = "Yes" if bool(poly_result) == bool(formula_result) else "No"
                    if match_str == "No":
                        all_match = False
                        
                    print(f"{int(a)} {int(b)} {int(c)} {int(d)} |      {int(poly_result)}      |       {int(formula_result)}        |  {match_str}")

    print("-" * 60)
    if all_match:
        print("\nVerification successful! The formula is equivalent to the polynomial.\n")
    else:
        print("\nVerification failed! The formula is not equivalent.\n")

    print("The original Zhigalkin polynomial with explicit coefficients is:")
    # The problem asks to "output each number in the final equation".
    # Since all coefficients in the given polynomial are 1, we print them explicitly.
    poly_str = "1*(d) ⊕ 1*(c) ⊕ 1*(b ∧ c) ⊕ 1*(a ∧ d) ⊕ 1*(a ∧ c) ⊕ 1*(a ∧ b ∧ d) ⊕ 1*(a ∧ b ∧ c)"
    print(poly_str)
    
    print("\nThe derived Boolean formula using the allowed operators is:")
    formula_str = "((a ∨ b) ∨ (c ↔ d)) → ¬(b → ¬d)"
    print(formula_str)

if __name__ == '__main__':
    solve()