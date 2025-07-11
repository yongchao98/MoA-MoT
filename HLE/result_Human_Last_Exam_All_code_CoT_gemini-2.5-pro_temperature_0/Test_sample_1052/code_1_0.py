import itertools

class ANF:
    """
    A class to represent and manipulate Algebraic Normal Form (ANF) polynomials.
    A polynomial is stored as a set of frozensets, where each frozenset
    represents a term (a product of variables).
    e.g., a + b*c + 1 is represented as {frozenset({'a'}), frozenset({'b', 'c'}), frozenset()}.
    """

    def __init__(self, terms=None):
        if terms is None:
            self.terms = set()
        else:
            self.terms = terms

    def __str__(self):
        if not self.terms:
            return "0"
        
        # Helper to format a single term
        def format_term(term):
            if not term:
                return "1"
            return "*".join(sorted(list(term)))

        # Sort terms for consistent output: by length, then alphabetically
        sorted_terms = sorted(list(self.terms), key=lambda t: (len(t), sorted(list(t))))
        
        return " + ".join(format_term(t) for t in sorted_terms)

    def __repr__(self):
        return f"ANF({self.terms})"

    def __eq__(self, other):
        return self.terms == other.terms

    def __add__(self, other):
        # Addition in GF(2) is XOR (symmetric difference of term sets)
        return ANF(self.terms.symmetric_difference(other.terms))

    def __xor__(self, other):
        # XOR is the same as addition
        return self.__add__(other)

    def __mul__(self, other):
        # Multiplication in GF(2) is AND
        new_terms = set()
        for t1 in self.terms:
            for t2 in other.terms:
                # The new term is the union of variables
                # The result is the XOR sum of all these products
                new_terms.symmetric_difference_update({t1.union(t2)})
        return ANF(new_terms)

# --- Helper functions to build formulas ---

def Var(name):
    """Creates a variable."""
    return ANF({frozenset({name})})

def Const(val):
    """Creates a constant (0 or 1)."""
    if val == 1:
        return ANF({frozenset()})
    return ANF()

# --- Logical operators converted to ANF operations ---

def anf_not(p):
    return Const(1) + p

def anf_or(p1, p2):
    return p1 + p2 + (p1 * p2)

def anf_and(p1, p2):
    return p1 * p2

def anf_implies(p1, p2):
    return Const(1) + p1 + (p1 * p2)

def anf_iff(p1, p2):
    return Const(1) + p1 + p2

def anf_nand(p1, p2):
    return Const(1) + (p1 * p2)

def anf_nor(p1, p2):
    return Const(1) + p1 + p2 + (p1 * p2)

def solve():
    """
    Solves the two-part problem.
    """
    print("--- Part 1: Verification ---")
    print("Goal: Verify that the ANF of (a -> b) is (((a & b) ^ a) ^ 1)\n")

    # Define variables
    a = Var('a')
    b = Var('b')

    # 1. Calculate ANF of (a -> b)
    formula_implies = anf_implies(a, b)
    print(f"The ANF of (a -> b) is derived as:")
    print(f"(a -> b) = (¬a ∨ b) = (1+a) ∨ b = (1+a) + b + (1+a)*b = 1+a+b+b+a*b = 1 + a + a*b")
    print(f"Calculated ANF of (a -> b): {formula_implies}\n")

    # 2. Simplify the given polynomial
    poly1 = (anf_and(a, b) ^ a) ^ Const(1)
    print(f"The given polynomial is (((a & b) ^ a) ^ 1)")
    print(f"Simplified polynomial: {poly1}\n")

    # 3. Compare
    if formula_implies == poly1:
        print("Verification successful: The two forms are equivalent.\n")
    else:
        print("Verification failed.\n")

    print("--- Part 2: Finding the Boolean Formula ---")
    print("Goal: Find the formula for the polynomial P2\n")
    
    # Define more variables
    c = Var('c')
    d = Var('d')

    # The given polynomial P2
    # P2 = ((((d ^ c) ^ (b & c)) ^ (a & d)) ^ (a & c)) ^ ((a & (b & d)) ^ (a & (b & c)))
    # Simplified: c ^ d ^ (b&c) ^ (a&d) ^ (a&c) ^ (a&b&d) ^ (a&b&c)
    p2_terms = {
        frozenset({'c'}), frozenset({'d'}), frozenset({'b', 'c'}), 
        frozenset({'a', 'd'}), frozenset({'a', 'c'}), 
        frozenset({'a', 'b', 'd'}), frozenset({'a', 'b', 'c'})
    }
    P2 = ANF(p2_terms)
    print(f"Target polynomial P2: {P2}\n")

    # Proposed candidate formula: ¬((c → (a ∨ b)) ↔ (d → ¬(a → b)))
    # Let's compute its ANF step-by-step
    print("Let's test the candidate formula: ¬((c → (a ∨ b)) ↔ (d → ¬(a → b)))")
    
    # Part A: c → (a ∨ b)
    a_or_b = anf_or(a, b)
    # print(f"ANF of (a ∨ b): {a_or_b}")
    part_A = anf_implies(c, a_or_b)
    print(f"ANF of (c → (a ∨ b)): {part_A}")

    # Part B: d → ¬(a → b)
    a_implies_b = anf_implies(a, b)
    not_a_implies_b = anf_not(a_implies_b)
    # print(f"ANF of ¬(a → b): {not_a_implies_b}")
    part_B = anf_implies(d, not_a_implies_b)
    print(f"ANF of (d → ¬(a → b)): {part_B}")

    # Combine with IFF (↔︎)
    combined_iff = anf_iff(part_A, part_B)
    print(f"ANF of ((c → (a ∨ b)) ↔ (d → ¬(a → b))): {combined_iff}")

    # Final negation (¬)
    final_formula_anf = anf_not(combined_iff)
    print(f"ANF of ¬((c → (a ∨ b)) ↔ (d → ¬(a → b))): {final_formula_anf}\n")

    # Check if it matches P2
    if final_formula_anf == P2:
        print("Success! The ANF of the candidate formula matches the target polynomial P2.")
        final_formula_str = "¬((c → (a ∨ b)) ↔ (d → ¬(a → b)))"
        print(f"\nThe derived Boolean formula is: {final_formula_str}")
        # This is the final answer for the prompt
        print(f"\n<<<{final_formula_str}>>>")
    else:
        print("Failure. The ANF of the candidate formula does not match P2.")

if __name__ == '__main__':
    solve()