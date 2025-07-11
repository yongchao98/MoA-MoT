def print_solution():
    """
    This function prints the solution to the problem.
    It first displays the full Zhigalkin polynomial equation.
    Then, it prints a derived Boolean formula that corresponds to this polynomial
    and uses the allowed set of logical operators.
    """
    
    # The given Zhigalkin polynomial
    zhigalkin_poly = "((((d \u2295 c) \u2295 (b \u2227 c)) \u2295 (a \u2227 d)) \u2295 (a \u2227 c)) \u2295 ((a \u2227 (b \u2227 d)) \u2295 (a \u2227 (b \u2227 c)))"
    
    # The derived Boolean formula using the allowed operators
    # This formula is logically equivalent to the one derived from the truth table:
    # (b and d) or ((not a and not b) and (c xnor d))
    # It has been converted to use the allowed operators.
    # Note: (not(b and d)) is equivalent to (b nand d), which is (b \u2191 d)
    # The structure is of the form X or Y, which can be written as (X -> Y) -> Y
    # Or more directly, not(X) -> Y
    
    # Let X = (b and d) and Y = ((a nor b) and (c xnor d))
    # The formula is X or Y, which can be expressed as not(X) -> Y
    # where not(X) = not(b and d) = b nand d = (b \u2191 d)
    # and Y = (a nor b) and (c xnor d) = (a \u2193 b) and (c \u2194 not d)
    # Y itself can be written as not((a \u2193 b) -> not(c \u2194 (not d)))
    boolean_formula = "\u00ac(\u00acb \u2228 \u00acd) \u2228 \u00ac(\u00ac(a \u2193 b) \u2228 \u00ac(c \u2194 \u00acd))"
    # A slightly more readable version using implication:
    # not(b and d) -> ((a nor b) and (c xnor d))
    # which expands to: (b \u2191 d) \u2192 \u00ac((a \u2193 b) \u2191 (c \u2194 \u00acd))
    
    # Let's use a version that is more directly translatable from the DNF.
    # F = (b \u2227 d) \u2228 (\u00ac a \u2227 \u00ac b \u2227 (c \u2194 \u00ac d))
    # \u00ac a \u2227 \u00ac b is equivalent to (a \u2193 b)
    # So F = (b \u2227 d) \u2228 ((a \u2193 b) \u2227 (c \u2194 \u00ac d))
    # This is a disjunction (OR) of two conjunctions (ANDs).
    # We can write X \u2227 Y as \u00ac(X \u2192 \u00acY)
    # Let's provide a final answer built with OR, NOT, NOR, and IFF, which are in the allowed list.
    final_formula = "\u00ac(\u00acb \u2228 \u00acd) \u2228 \u00ac(\u00ac(a \u2193 b) \u2228 \u00ac(c \u2194 \u00acd))"

    print("The Zhigalkin polynomial is:")
    print(zhigalkin_poly)
    print("\nA Boolean formula from which the polynomial is derived is:")
    print(final_formula)

print_solution()
<<<(b ∧ d) ∨ (¬a ∧ ¬b ∧ (c ↔ ¬d))>>>