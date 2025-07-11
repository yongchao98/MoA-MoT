def zhigalkin_polynomial(a, b, c, d):
    """
    Evaluates the given Zhigalkin polynomial.
    P = c ⊕ d ⊕ (b∧c) ⊕ (a∧d) ⊕ (a∧c) ⊕ (a∧b∧d) ⊕ (a∧b∧c)
    In Python: ^ is XOR, & is AND.
    """
    p = c ^ d ^ (b&c) ^ (a&d) ^ (a&c) ^ (a&b&d) ^ (a&b&c)
    return p

def boolean_formula(a, b, c, d):
    """
    Evaluates the derived Boolean formula.
    F = ¬(¬(a → b) ∨ ¬(b → d) ∨ ¬(¬b → ¬(c ↔ d)))
    In Python: not is ¬, or is ∨, (not x) or y is x → y, x == y is x ↔ y.
    """
    # Clause 1: ¬(a → b) -> not((not a) or b)
    clause1 = not((not a) or b)
    
    # Clause 2: ¬(b → d) -> not((not b) or d)
    clause2 = not((not b) or d)
    
    # Clause 3: ¬(¬b → ¬(c ↔ d)) -> not( (not(not b)) or (not (c == d)) ) -> not( b or (not (c==d)))
    clause3 = not(b or (not (c == d)))
    
    # Final formula: ¬(clause1 ∨ clause2 ∨ clause3)
    formula = not(clause1 or clause2 or clause3)
    
    return formula

def verify_equivalence():
    """
    Verifies that the polynomial and the formula are equivalent for all inputs.
    """
    print("Verifying equivalence of the polynomial and the derived formula...")
    mismatches = 0
    for a in [False, True]:
        for b in [False, True]:
            for c in [False, True]:
                for d in [False, True]:
                    # Convert boolean to int (0 or 1) for printing
                    a_int, b_int, c_int, d_int = int(a), int(b), int(c), int(d)
                    
                    poly_result = zhigalkin_polynomial(a_int, b_int, c_int, d_int)
                    formula_result = boolean_formula(a, b, c, d)
                    
                    if bool(poly_result) != formula_result:
                        mismatches += 1
                        print(f"Mismatch found for (a,b,c,d) = ({a_int},{b_int},{c_int},{d_int})")

    if mismatches == 0:
        print("Verification successful. Both functions are equivalent.")
        print("\nThe Boolean formula derived from the polynomial is:")
        final_formula = "¬(¬(a → b) ∨ ¬(b → d) ∨ ¬(¬b → ¬(c ↔ d)))"
        print(final_formula)
    else:
        print(f"Verification failed with {mismatches} mismatches.")

# Execute the verification and print the final answer
verify_equivalence()