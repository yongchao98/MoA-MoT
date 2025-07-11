def zhigalkin_polynomial(a, b, c, d):
    """
    Evaluates the given Zhigalkin polynomial for boolean inputs.
    Polynomial: ((((d^c)^(b&c))^(a&d)^(a&c))^((a&b&d)^(a&b&c)))
    Note: In Python, a^b is XOR, a&b is AND. Boolean inputs are treated as 0/1.
    """
    term1 = d ^ c
    term2 = b & c
    term3 = a & d
    term4 = a & c
    term5 = a & b & d
    term6 = a & b & c
    
    # XOR is associative, so we can XOR all terms together
    result = term1 ^ term2 ^ term3 ^ term4 ^ term5 ^ term6
    return result

def boolean_formula(a, b, c, d):
    """
    Evaluates the derived Boolean formula.
    Formula: ((a ∨ b) ∨ (c ↔ d)) → ¬(b ↑ d)
    Note: Python's 'or', '==', 'not', and '(not (a and b))' correspond to ∨, ↔, ¬, and ↑.
    """
    # Evaluate sub-expressions
    b_nand_d = not (b and d)  # b ↑ d
    not_b_nand_d = not b_nand_d # ¬(b ↑ d)
    
    c_xnor_d = (c == d)  # c ↔ d
    a_or_b = a or b      # a ∨ b
    
    left_side = a_or_b or c_xnor_d # (a ∨ b) ∨ (c ↔ d)
    
    # Implication A → B is equivalent to (¬A ∨ B)
    result = (not left_side) or not_b_nand_d
    return result

def solve():
    """
    Verifies the equivalence of the polynomial and the formula and prints the results.
    """
    print("Verifying that the Zhigalkin polynomial and the Boolean formula are equivalent.")
    print("a b c d | Polynomial | Formula | Match?")
    print("---------------------------------------")
    
    all_match = True
    for a_val in [False, True]:
        for b_val in [False, True]:
            for c_val in [False, True]:
                for d_val in [False, True]:
                    poly_res = zhigalkin_polynomial(a_val, b_val, c_val, d_val)
                    form_res = boolean_formula(a_val, b_val, c_val, d_val)
                    
                    match_status = "Yes" if poly_res == form_res else "No"
                    if not (poly_res == form_res):
                        all_match = False

                    print(f"{int(a_val)} {int(b_val)} {int(c_val)} {int(d_val)} |      {int(poly_res)}     |    {int(form_res)}    |   {match_status}")

    if all_match:
        print("\nVerification successful! All 16 outputs match.")
    else:
        print("\nVerification failed! Mismatch found.")

    final_formula = "((a ∨ b) ∨ (c ↔ d)) → ¬(b ↑ d)"
    print("\nThe derived Boolean formula is:")
    print(final_formula)

    print("\nThe original polynomial equation implicitly contains numbers (coefficients for each term).")
    print("Polynomial: (1*d ⊕ 1*c) ⊕ 1*(b∧c) ⊕ 1*(a∧d) ⊕ 1*(a∧c) ⊕ 1*(a∧b∧d) ⊕ 1*(a∧b∧c)")
    print("The numbers (coefficients) in this equation are:")
    # There are 7 terms being XORed together, each with an implicit coefficient of 1.
    coeffs = [1, 1, 1, 1, 1, 1, 1]
    for num in coeffs:
        print(num)

if __name__ == "__main__":
    solve()