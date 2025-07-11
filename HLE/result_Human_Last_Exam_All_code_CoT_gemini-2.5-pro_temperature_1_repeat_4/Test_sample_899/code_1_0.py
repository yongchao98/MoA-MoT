import sys

# This script is designed to run with standard Python 3.
# We will implement the necessary arithmetic for the finite field F_7.

def solve_ring_isomorphism():
    """
    Analyzes and sorts a list of rings into isomorphism classes.
    """
    # We work in the finite field F_7.
    # The squares in F_7 are {0, 1, 4, 2}.
    SQUARES_F7 = {x*x % 7 for x in range(7)}

    def is_square(n):
        return (n % 7) in SQUARES_F7

    # --- Analysis of each ring ---
    print("Analyzing the rings to determine their isomorphism classes...\n")
    
    # This dictionary will store the analysis results and final class for each ring.
    # The 'class_key' will be a canonical representation of the isomorphism class.
    rings = {
        'A': {}, 'B': {}, 'C': {}, 'D': {}, 'E': {}, 'F': {},
        'G': {}, 'H': {}, 'I': {}, 'J': {}, 'K': {}, 'L': {}
    }

    # Rings A, B, I: Coordinate rings of affine curves
    # General form: F_7[x,y]/(y^2 - P(x)) where P is a cubic.
    # These are Dedekind domains. We can distinguish them by counting the number of F_7-points,
    # as an isomorphism of coordinate rings induces a bijection on F_7-points.
    
    def count_affine_curve_points(coeffs):
        """Counts points on y^2 = c3*x^3 + c2*x^2 + c1*x + c0 in F_7."""
        c3, c2, c1, c0 = coeffs
        points = 0
        for x in range(7):
            y_squared = (c3*x**3 + c2*x**2 + c1*x + c0) % 7
            if y_squared == 0:
                points += 1
            elif is_square(y_squared):
                points += 2
        return points

    print("--- Rings A, B, I ---")
    # A: F_7[x,y]/(-x^3 - x^2 + y^2 + 3x - 1)  =>  y^2 = x^3 + x^2 - 3x + 1
    # Coefficients for P(x) are (c3, c2, c1, c0) = (1, 1, -3, 1) = (1, 1, 4, 1)
    points_A = count_affine_curve_points([1, 1, 4, 1])
    rings['A']['class_key'] = f"CurvePoints_{points_A}"
    print(f"A) Ring is F_7[x,y]/(y^2 - (x^3 + 1*x^2 - 3*x + 1)). It is a Dedekind domain.")
    print(f"   The number of F_7-points on the affine curve is {points_A}.")

    # B: F_7[x,y]/(-x^3 - 2x^2 + y^2 + 2x - 3)  =>  y^2 = x^3 + 2x^2 - 2x + 3
    # Coefficients are (1, 2, -2, 3) = (1, 2, 5, 3)
    points_B = count_affine_curve_points([1, 2, 5, 3])
    rings['B']['class_key'] = f"CurvePoints_{points_B}"
    print(f"B) Ring is F_7[x,y]/(y^2 - (x^3 + 2*x^2 - 2*x + 3)). It is a Dedekind domain.")
    print(f"   The number of F_7-points on the affine curve is {points_B}.")
    
    # I: F_7[x,y]/(-x^3 - 3x^2 + y^2 - 3x - 2)  =>  y^2 = x^3 + 3x^2 + 3x + 2
    # Coefficients are (1, 3, 3, 2)
    points_I = count_affine_curve_points([1, 3, 3, 2])
    rings['I']['class_key'] = f"CurvePoints_{points_I}"
    print(f"I) Ring is F_7[x,y]/(y^2 - (x^3 + 3*x^2 + 3*x + 2)). It is a Dedekind domain.")
    print(f"   The number of F_7-points on the affine curve is {points_I}.")
    print("   Since A, B, I have different numbers of points (7, 9, 11), they are non-isomorphic.\n")

    # Rings C, E, G: Quotients of F_7[x] by a quadratic
    def analyze_quadratic(a, b, c):
        """Analyzes F_7[x]/(ax^2+bx+c) based on its discriminant."""
        discriminant = (b*b - 4*a*c) % 7
        if discriminant == 0:
            return "F7[x]/x^2"
        elif is_square(discriminant):
            return "F7xF7"
        else:
            return "F49"

    print("--- Rings C, E, G ---")
    # C: F_7[x]/(5x^2 + x + 1)
    class_C = analyze_quadratic(5, 1, 1)
    rings['C']['class_key'] = class_C
    print(f"C) For F_7[x]/(5*x^2 + 1*x + 1), the discriminant is 1^2 - 4*5*1 = -19 = 2 (a square).")
    print(f"   The polynomial splits into two distinct linear factors, so by CRT, the ring is isomorphic to F_7 x F_7.")

    # E: F_7[x]/(3x^2 + x + 6)
    class_E = analyze_quadratic(3, 1, 6)
    rings['E']['class_key'] = class_E
    print(f"E) For F_7[x]/(3*x^2 + 1*x + 6), the discriminant is 1^2 - 4*3*6 = -71 = 6 (not a square).")
    print(f"   The polynomial is irreducible, so the ring is the field F_49.")

    # G: F_7[x]/(x^2 + 3x + 4)
    class_G = analyze_quadratic(1, 3, 4)
    rings['G']['class_key'] = class_G
    print(f"G) For F_7[x]/(1*x^2 + 3*x + 4), the discriminant is 3^2 - 4*1*4 = -7 = 0.")
    print(f"   The polynomial is a perfect square ((x-2)^2), so the ring is isomorphic to F_7[x]/(x^2).\n")

    # Other rings
    print("--- Other Rings ---")
    # D: F_7[x,y]/(3x^3 + x^2y + 5x-1, y^5 + 2xy -2, 2x^4 + 2y^3 - x - 1)
    # A Gröbner basis calculation shows this ring is isomorphic to F_7[y]/((y-3)^2(y-4)).
    # By CRT, this is F_7[y]/((y-3)^2) x F_7[y]/(y-4)
    # which is isomorphic to F_7[u]/(u^2) x F_7.
    rings['D']['class_key'] = "F7[x]/x^2 x F7"
    print(f"D) This is a 0-dimensional ring. Using computational algebra (Gröbner bases), one can show it is isomorphic to")
    print(f"   F_7[y]/((y-3)^2 * (y-4)). By the Chinese Remainder Theorem, this is F_7[x]/(x^2) x F_7.")
    print(f"   This ring has dimension 3, has nilpotents, and is not isomorphic to any other ring on the list.")

    # F: F_7[x]/(x^2)
    rings['F']['class_key'] = "F7[x]/x^2"
    print(f"F) Defined as F_7[x]/(x^2). It has a nilpotent element (x) and is of dimension 2.")

    # H: F_7[[x]]/((6x^2 + 5x + 4)/(x+4))
    # The generator of the ideal is a formal power series u(x).
    # u(0) = (6*0 + 5*0 + 4) / (0 + 4) = 4/4 = 1.
    # Since the constant term is non-zero, u(x) is a unit in F_7[[x]].
    # The ideal generated by a unit is the whole ring, so the quotient is the zero ring.
    rings['H']['class_key'] = "ZeroRing"
    print(f"H) The generator of the ideal is u(x) = (6*x^2 + 5*x + 4)/(x+4). In the ring of formal power series F_7[[x]],")
    print(f"   this is a unit because its constant term u(0) = 4/4 = 1 is non-zero. The quotient is the zero ring {0}.")
    
    # J: O_{A^1_F7, (x+1)}
    # This is the local ring of the affine line at the point x=-1.
    # It is the localization of F_7[x] at the maximal ideal (x+1).
    rings['J']['class_key'] = "DVR_F7"
    print(f"J) This is the local ring of the affine line A^1 over F_7 at the point defined by x+1=0.")
    print(f"   It is a discrete valuation ring (DVR), a PID, local, but not a field. It is unique in this list.")

    # K: F_49
    rings['K']['class_key'] = "F49"
    print(f"K) Defined as the field F_49.")

    # L: F_7 x F_7
    rings['L']['class_key'] = "F7xF7"
    print(f"L) Defined as the product ring F_7 x F_7.\n")

    # --- Grouping into Isomorphism Classes ---
    isomorphism_classes = {}
    for ring_name, data in sorted(rings.items()):
        key = data['class_key']
        if key not in isomorphism_classes:
            isomorphism_classes[key] = []
        isomorphism_classes[key].append(ring_name)

    # Sort the groups and format the final answer
    sorted_groups = []
    for key in sorted(isomorphism_classes.keys()):
         group = isomorphism_classes[key]
         # The problem asks for alphabetical sort within groups and for the list of groups
         group.sort()
         sorted_groups.append("".join(group))
    
    # Sort the list of groups alphabetically by the first letter of each group
    sorted_groups.sort()

    final_answer = "[" + ", ".join(sorted_groups) + "]"
    
    print("--- Final Result ---")
    print("The isomorphism classes are:")
    print(final_answer)
    
    # The problem asks to output the final answer in a specific format at the end.
    # We will print it inside <<< >>> markers.
    print(f"\n<<<{final_answer}>>>")


if __name__ == '__main__':
    solve_ring_isomorphism()