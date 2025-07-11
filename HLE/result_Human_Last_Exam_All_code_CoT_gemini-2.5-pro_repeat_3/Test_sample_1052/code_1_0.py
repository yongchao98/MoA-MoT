import itertools

def b_and(x, y):
    """Boolean AND"""
    return x & y

def b_xor(x, y):
    """Boolean XOR"""
    return x ^ y

def b_or(x, y):
    """Boolean OR"""
    return x | y

def b_implies(x, y):
    """Boolean Implication (x -> y)"""
    return (1 - x) | y

def b_iff(x, y):
    """Boolean Biconditional / IFF (x <-> y)"""
    return 1 if x == y else 0

def evaluate_polynomial(a, b, c, d):
    """
    Evaluates the Zhigalkin polynomial:
    ((((d ⊕ c) ⊕ (b ∧ c)) ⊕ (a ∧ d)) ⊕ (a ∧ c)) ⊕ ((a ∧ (b ∧ d)) ⊕ (a ∧ (b ∧ c)))
    which simplifies to:
    c ⊕ d ⊕ (b∧c) ⊕ (a∧d) ⊕ (a∧c) ⊕ (a∧b∧d) ⊕ (a∧b∧c)
    """
    # Grouping terms for calculation
    p = b_xor(d, c)
    p = b_xor(p, b_and(b, c))
    p = b_xor(p, b_and(a, d))
    p = b_xor(p, b_and(a, c))
    p = b_xor(p, b_and(a, b_and(b, d)))
    p = b_xor(p, b_and(a, b_and(b, c)))
    return p

def evaluate_formula(a, b, c, d):
    """
    Evaluates the derived Boolean formula:
    (a ∨ b ∨ (c ↔ d)) → (b ∧ d)
    """
    # Left-hand side of implication
    lhs_iff = b_iff(c, d)
    lhs_or = b_or(b_or(a, b), lhs_iff)

    # Right-hand side of implication
    rhs_and = b_and(b, d)
    
    return b_implies(lhs_or, rhs_and)

# --- Main Verification Logic ---

# The polynomial and the derived formula strings
polynomial_str = "((((d ⊕ c) ⊕ (b ∧ c)) ⊕ (a ∧ d)) ⊕ (a ∧ c)) ⊕ ((a ∧ (b ∧ d)) ⊕ (a ∧ (b ∧ c)))"
derived_formula_str = "(a ∨ b ∨ (c ↔ d)) → (b ∧ d)"

print("This script verifies that the given Zhigalkin polynomial is equivalent to the derived Boolean formula.")
print("-" * 70)
print(f"Zhigalkin Polynomial P(a,b,c,d):\n  {polynomial_str}\n")
print(f"Derived Boolean Formula F(a,b,c,d):\n  {derived_formula_str}")
print("-" * 70)
print("Verification by comparing truth tables for all 16 inputs:")
print("\n  a b c d | P | F | Match?")
print("  --------------------------")

all_match = True
# Iterate through all 2^4 = 16 combinations of inputs
for a, b, c, d in itertools.product([0, 1], repeat=4):
    poly_result = evaluate_polynomial(a, b, c, d)
    formula_result = evaluate_formula(a, b, c, d)
    
    match_status = "✅" if poly_result == formula_result else "❌"
    if poly_result != formula_result:
        all_match = False
        
    print(f"  {a} {b} {c} {d} | {poly_result} | {formula_result} |   {match_status}")

print("  --------------------------")

if all_match:
    print("\nConclusion: The verification was successful. The derived formula is correct.")
else:
    print("\nConclusion: Verification failed. The formula is not equivalent to the polynomial.")

print("\nThe final Boolean formula derived from the polynomial is:")
# The final equation contains variables and operators, which are printed here.
print(derived_formula_str)
