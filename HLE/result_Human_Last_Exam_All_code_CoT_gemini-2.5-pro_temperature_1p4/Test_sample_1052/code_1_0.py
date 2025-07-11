def target_polynomial(a, b, c, d):
    """Calculates the value of the given Zhigalkin polynomial."""
    # Using integer arithmetic where XOR is addition mod 2, and AND is multiplication.
    # Polynomial: c ⊕ d ⊕ (b∧c) ⊕ (a∧d) ⊕ (a∧c) ⊕ (a∧b∧d) ⊕ (a∧b∧c)
    val = (c + d + (b*c) + (a*d) + (a*c) + (a*b*d) + (a*b*c)) % 2
    return val

def final_formula(a, b, c, d):
    """
    Calculates the value of the derived Boolean formula.
    Formula: (a ↑ (b → ¬d)) ↔ (a ∨ ((c → b) ↔ d))
    """
    # Helper functions for logical operations on integers (0/1)
    op_not = lambda x: 1 - x
    op_or = lambda x, y: min(x + y, 1)
    op_imp = lambda x, y: op_or(op_not(x), y)
    op_iff = lambda x, y: 1 if x == y else 0
    op_nand = lambda x, y: op_not(x * y)

    # Evaluate the right-hand side of the main equivalence
    # a ∨ ((c → b) ↔ d)
    c_imp_b = op_imp(c, b)
    cib_iff_d = op_iff(c_imp_b, d)
    rhs = op_or(a, cib_iff_d)
    
    # Evaluate the left-hand side of the main equivalence
    # a ↑ (b → ¬d)
    not_d = op_not(d)
    b_imp_not_d = op_imp(b, not_d)
    lhs = op_nand(a, b_imp_not_d)

    # Evaluate the main equivalence
    result = op_iff(lhs, rhs)
    return result

# Verification process
is_correct = True
for a_val in [0, 1]:
    for b_val in [0, 1]:
        for c_val in [0, 1]:
            for d_val in [0, 1]:
                poly_res = target_polynomial(a_val, b_val, c_val, d_val)
                formula_res = final_formula(a_val, b_val, c_val, d_val)
                if poly_res != formula_res:
                    print(f"Error: Mismatch at (a={a_val}, b={b_val}, c={c_val}, d={d_val})")
                    is_correct = False
                    break
            if not is_correct: break
        if not is_correct: break
    if not is_correct: break

if is_correct:
    # The derived formula is correct. Print it as the final answer.
    # The operators are: ¬ (NOT), ↑ (NAND), ↔ (IFF/XNOR), → (IMPLIES), ∨ (OR)
    final_equation = "(a ↑ (b → ¬d)) ↔ (a ∨ ((c → b) ↔ d))"
    print("The derived Boolean formula is:")
    print(final_equation)

<<<(a ↑ (b → ¬d)) ↔ (a ∨ ((c → b) ↔ d))>>>