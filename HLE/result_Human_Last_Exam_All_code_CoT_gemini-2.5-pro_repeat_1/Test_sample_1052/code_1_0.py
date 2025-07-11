def original_polynomial(a, b, c, d):
    """
    This function evaluates the given Zhigalkin polynomial.
    In Python, XOR is '^' and AND is '&'.
    """
    # The polynomial is:
    # c ⊕ d ⊕ (b∧c) ⊕ (a∧d) ⊕ (a∧c) ⊕ (a∧b∧d) ⊕ (a∧b∧c)
    p = (
        c ^ d ^
        (b & c) ^
        (a & d) ^
        (a & c) ^
        (a & b & d) ^
        (a & b & c)
    )
    return p

def derived_formula(a, b, c, d):
    """
    This function evaluates the Boolean formula derived from analyzing the polynomial's truth table.
    The formula is (¬b ∨ ¬d) → ¬(¬(a ∨ b) → (c ↔ d)).
    We translate it to Python's logic operators.
    a → b  is equivalent to (not a) or b
    a ↔ b  is equivalent to a == b
    """
    # Let's break down the formula: (¬b ∨ ¬d) → ¬(¬(a ∨ b) → (c ↔ d))
    
    # Part 1: (¬b ∨ ¬d)
    antecedent = (not b) or (not d)
    
    # Part 2: ¬(¬(a ∨ b) → (c ↔ d))
    # Part 2a: ¬(a ∨ b)
    inner_antecedent = not (a or b)
    # Part 2b: c ↔ d
    inner_consequent = (c == d)
    # Part 2c: ¬(a ∨ b) → (c ↔ d)
    inner_implication = (not inner_antecedent) or inner_consequent
    # Part 2d: ¬(...)
    consequent = not inner_implication
    
    # Final implication
    result = (not antecedent) or consequent
    return int(result)

def solve():
    """
    Verifies that the polynomial and the derived formula are equivalent
    and prints the final formula.
    """
    is_equivalent = True
    for a in [0, 1]:
        for b in [0, 1]:
            for c in [0, 1]:
                for d in [0, 1]:
                    if original_polynomial(a, b, c, d) != derived_formula(a, b, c, d):
                        is_equivalent = False
                        print(f"Mismatch found for a={a}, b={b}, c={c}, d={d}")
                        break
                if not is_equivalent:
                    break
            if not is_equivalent:
                break
        if not is_equivalent:
            break

    if is_equivalent:
        # The problem statement mentions "output each number in the final equation".
        # The given polynomial can be seen as a sum of terms with coefficient 1.
        # P = 1*(d⊕c) ⊕ 1*(b∧c) ⊕ 1*(a∧d) ⊕ 1*(a∧c) ⊕ 1*(a∧b∧d) ⊕ 1*(a∧b∧c)
        # There are no other numbers, so we will just print the final formula.
        final_formula = "(¬b ∨ ¬d) → ¬(¬(a ∨ b) → (c ↔ d))"
        print("The Boolean formula derived from the polynomial is:")
        print(final_formula)
        
        # Output the final answer in the requested format
        print(f"\n<<<{final_formula}>>>")

solve()