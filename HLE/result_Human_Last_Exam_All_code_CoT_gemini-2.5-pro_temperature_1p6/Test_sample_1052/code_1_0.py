# Helper function to perform XOR sum on dictionary-based polynomials
def poly_xor(p1, p2):
    res = p1.copy()
    for term, coeff in p2.items():
        res[term] = res.get(term, 0) ^ coeff
    # Remove terms with zero coefficient
    return {t: c for t, c in res.items() if c != 0}

# Helper function to perform AND (multiplication) on polynomials
def poly_mul(p1, p2):
    res = {}
    for t1, c1 in p1.items():
        for t2, c2 in p2.items():
            # Combine terms by sorting characters to handle commutativity (e.g., 'ab' == 'ba')
            term = "".join(sorted(set(t1 + t2)))
            # In GF(2), coefficients are just XORed
            res[term] = res.get(term, 0) ^ (c1 * c2)
    return {t: c for t, c in res.items() if c != 0}

# Helper to format and print a polynomial
def print_poly(name, p):
    if not p:
        print(f"P({name}) = 0")
        return
    # Sort terms for consistent output: by length, then alphabetically
    sorted_terms = sorted(p.keys(), key=lambda t: (len(t), t))
    poly_str = " \u2295 ".join(f"{p[t]}*{t}" if t else str(p[t]) for t in sorted_terms if p[t] != 0).replace("1*", "")
    print(f"P({name}) = {poly_str}")

# --- Main logic ---
# Base polynomials for variables
p_a = {'a': 1}
p_b = {'b': 1}
p_c = {'c': 1}
p_d = {'d': 1}
p_1 = {'': 1} # Constant 1 for True

print("Let's derive the Zhigalkin polynomial for the formula F = not((c and (a nor b))leftrightarrow(d and (a -> b)))")
print("This will show that it matches the polynomial given in the question.\n")

# Step 1: Calculate P(a -> b)
# a -> b == not(a) or b ==> poly is 1 + a + ab
print("Step 1: Calculate the polynomial for (a -> b)")
p_a_implies_b = poly_xor(poly_xor(p_1, p_a), poly_mul(p_a, p_b))
print_poly("a -> b", p_a_implies_b)
print("-" * 30)

# Step 2: Calculate P(a nor b)
# a nor b == not(a or b) ==> poly is 1 + a + b + ab
print("Step 2: Calculate the polynomial for (a nor b)")
p_a_nor_b = poly_xor(p_a_implies_b, p_b) # Quick way since (a->b) = (a nor b) xor b
print_poly("a nor b", p_a_nor_b)
print("-" * 30)

# Step 3: Calculate P(c and (a nor b))
print("Step 3: Calculate the polynomial for A = (c and (a nor b))")
p_A = poly_mul(p_c, p_a_nor_b)
print_poly("A", p_A)
print("-" * 30)

# Step 4: Calculate P(d and (a -> b))
print("Step 4: Calculate the polynomial for B = (d and (a -> b))")
p_B = poly_mul(p_d, p_a_implies_b)
print_poly("B", p_B)
print("-" * 30)

# Step 5: Calculate P(A <-> B)
# A <-> B ==> poly is 1 + P(A) + P(B)
print("Step 5: Calculate the polynomial for (A <-> B)")
p_A_xnor_B = poly_xor(poly_xor(p_1, p_A), p_B)
print_poly("A <-> B", p_A_xnor_B)
print("-" * 30)

# Step 6: Calculate P(not(A <-> B))
# not(Z) ==> poly is 1 + P(Z)
print("Step 6: Calculate the final polynomial for not(A <-> B)")
p_final = poly_xor(p_1, p_A_xnor_B)
print("The final polynomial is:")
print_poly("F", p_final)
print("\nThis polynomial c \u2295 d \u2295 ac \u2295 ad \u2295 bc \u2295 abc \u2295 abd is a reordered version of the one from the problem statement.")

# We don't need to replace `and` for the calculation, but the final formula must use the allowed ops.
final_formula_str = "not((not(c nand (a nor b)))leftrightarrow(not(d nand (a -> b))))"
print(f"\nA valid boolean formula is: {final_formula_str}")
