def original_polynomial(a, b, c, d):
    """
    Calculates the value of the given Zhigalkin polynomial.
    P = ((((d ^ c) ^ (b & c)) ^ (a & d)) ^ (a & c)) ^ ((a & (b & d)) ^ (a & (b & c)))
    """
    # XOR is associative, so we can remove the parentheses for a flat expression.
    # Note: In Python, `^` is XOR and `&` is AND.
    return d ^ c ^ (b & c) ^ (a & d) ^ (a & c) ^ (a & b & d) ^ (a & b & c)

def derived_formula(a, b, c, d):
    """
    Calculates the value of the derived Boolean formula.
    Formula: not(not(a -> b) or (c <-> d))
    """
    # Boolean equivalents for logical operators
    # a -> b is equivalent to (not a) or b
    # a <-> b is equivalent to (a == b) or not (a ^ b)
    
    val_a_implies_b = (not a) or b
    val_c_equiv_d = (a == a) if (c == d) else (a != a) # A simple way to get True/False
    
    # The formula is not(not(a -> b) or (c <-> d))
    return not (not val_a_implies_b or val_c_equiv_d)

def verify_and_print_solution():
    """
    Verifies that the polynomial and the formula are equivalent and prints the result.
    """
    all_match = True
    for a_val in [False, True]:
        for b_val in [False, True]:
            for c_val in [False, True]:
                for d_val in [False, True]:
                    # Convert boolean to int (0 or 1) for polynomial calculation
                    a_int, b_int, c_int, d_int = int(a_val), int(b_val), int(c_val), int(d_val)
                    
                    poly_result = bool(original_polynomial(a_int, b_int, c_int, d_int))
                    formula_result = derived_formula(a_val, b_val, c_val, d_val)
                    
                    if poly_result != formula_result:
                        print(f"Mismatch found for a={a_int}, b={b_int}, c={c_int}, d={d_int}")
                        all_match = False
                        break
                if not all_match: break
            if not all_match: break
        if not all_match: break

    if all_match:
        print("Verification successful: The polynomial and the derived formula are equivalent.")
        print("\nA Boolean formula for the given polynomial is:")
        # The prompt asks to output numbers in the equation, but there are no constants (like 0 or 1)
        # in the final simplified formula. We will print the formula structure itself.
        print("¬(¬(a → b) ∨ (c ↔ d))")

# Execute the verification and print the solution
verify_and_print_solution()
