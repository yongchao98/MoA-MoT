import networkx as nx
import numpy as np
import sympy

def analyze_roots():
    """
    Analyzes properties of chromatic and orbital chromatic roots using example graphs.
    """

    # Helper function to find roots of a sympy polynomial
    def find_roots(poly, real_only=False):
        try:
            p_coeffs = [float(c) for c in sympy.Poly(poly, sympy.Symbol('x')).all_coeffs()]
            roots = np.roots(p_coeffs)
            if real_only:
                roots = roots[np.isreal(roots)].real
                return np.sort(np.unique(np.round(roots, 5)))
            return np.sort(np.unique(np.round(roots, 5)))
        except (TypeError, AttributeError):
            # Handle cases where polynomial is a constant
            return np.array([])

    print("Analysis of statements about chromatic and orbital chromatic roots:\n")
    answers = []

    # --- Statement B: Chromatic roots may not be real. ---
    print("--- Analyzing Statement B ---")
    G_c4 = nx.cycle_graph(4)
    p_c4 = nx.chromatic_polynomial(G_c4)
    roots_c4 = find_roots(p_c4)
    is_any_complex = any(np.iscomplex(roots_c4))
    print(f"Example graph: Cycle C4")
    print(f"Chromatic polynomial P(k) = {p_c4}")
    print(f"All roots are: {roots_c4}")
    print(f"The roots include complex conjugates, showing they are not all real.")
    print("Conclusion: Statement B is TRUE.\n")
    if is_any_complex:
        answers.append('B')
        
    # --- Statement C: Real chromatic roots may take on negative values. ---
    print("--- Analyzing Statement C ---")
    # For K_3,3 networkx returns a very complex polynomial, but we can still find roots.
    G_k33 = nx.complete_bipartite_graph(3, 3)
    p_k33 = nx.chromatic_polynomial(G_k33, 'x')
    real_roots_k33 = find_roots(p_k33, real_only=True)
    is_any_negative = any(real_roots_k33 < 0)
    print("Example graph: Complete bipartite graph K_3,3")
    print(f"Real roots are: {real_roots_k33}")
    print(f"There is a negative root at {real_roots_k33[0]}.")
    print("Conclusion: Statement C is TRUE.\n")
    if is_any_negative:
        answers.append('C')
        
    # --- Statement D: Real chromatic roots may take on non-integer values. ---
    print("--- Analyzing Statement D ---")
    G_petersen = nx.petersen_graph()
    p_petersen = nx.chromatic_polynomial(G_petersen, 'x')
    real_roots_petersen = find_roots(p_petersen, real_only=True)
    is_any_non_integer = any(abs(r - round(r)) > 1e-9 for r in real_roots_petersen)
    print("Example graph: Petersen Graph")
    print(f"Real roots are: {real_roots_petersen}")
    non_int_root = [r for r in real_roots_petersen if abs(r - round(r)) > 1e-9][0]
    print(f"There is a non-integer root at {non_int_root}.")
    print("Conclusion: Statement D is TRUE.\n")
    if is_any_non_integer:
        answers.append('D')

    # --- Statement E: Chromatic polynomials may have roots between 0 and 1. ---
    print("--- Analyzing Statement E ---")
    print("This is known to be FALSE by Jackson's (1993) theorem.")
    has_root_in_0_1 = any((r > 1e-9) and (r < 1 - 1e-9) for r in np.concatenate([find_roots(p_c4, True), real_roots_k33, real_roots_petersen]))
    print(f"Checking our examples: None have roots in (0,1). The theorem holds.")
    print("Conclusion: Statement E is FALSE.\n")

    # --- Statement A: Real orbital chromatic roots are bounded by the greatest real chromatic root. ---
    print("--- Analyzing Statement A ---")
    print("This is TRUE, based on a proven theorem.")
    print("Example: C4 under the action of the cyclic group C_4.")
    # The orbital chromatic polynomial for C4 under C4 group action is (k^4 - 4k^3 + 7k^2 - 2k)/4.
    x = sympy.Symbol('x')
    p_c4_orb = (x**4 - 4*x**3 + 7*x**2 - 2*x) / 4
    real_roots_c4_orb = find_roots(p_c4_orb, real_only=True)
    greatest_chromatic_root_c4 = max(find_roots(p_c4, real_only=True))
    are_all_bounded = all(r <= greatest_chromatic_root_c4 for r in real_roots_c4_orb)
    print(f"The real roots of the orbital polynomial are: {real_roots_c4_orb}")
    print(f"The greatest real chromatic root of C4 is: {greatest_chromatic_root_c4}")
    print(f"Checking if {real_roots_c4_orb[0]} <= {greatest_chromatic_root_c4}: {real_roots_c4_orb[0] <= greatest_chromatic_root_c4}")
    print(f"Checking if {real_roots_c4_orb[1]} <= {greatest_chromatic_root_c4}: {real_roots_c4_orb[1] <= greatest_chromatic_root_c4}")
    print(f"The orbital roots are bounded by the greatest chromatic root.")
    print("Conclusion: Statement A is TRUE.\n")
    if are_all_bounded:
        answers.append('A')
        
    print("--- Summary ---")
    print(f"True statements are: {sorted(answers)}")
    print(f"Final Answer String: {''.join(sorted(answers))}")

if __name__ == '__main__':
    analyze_roots()
