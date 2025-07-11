import math

def v2(n):
    """
    Computes the 2-adic valuation of an integer n.
    Returns infinity if n is 0.
    """
    if n == 0:
        return float('inf')
    if n % 2 != 0:
        return 0
    return (n & -n).bit_length() - 1

def print_polynomial(coeffs_dict, var='x'):
    """Prints a polynomial from a dictionary of coefficients."""
    terms = []
    for power, coeff in sorted(coeffs_dict.items(), reverse=True):
        if coeff == 0:
            continue
        if power == 0:
            terms.append(f"{coeff}")
        elif power == 1:
            terms.append(f"{coeff}*{var}")
        else:
            terms.append(f"{coeff}*{var}^{power}")
    return " + ".join(terms)

def main():
    print("Step 1: Analyze the curve and its branch points.")
    p_coeffs = {5: 8, 4: 4, 3: 4, 2: 1, 1: 8, 0: 0}
    print(f"The curve is y^2 = {print_polynomial(p_coeffs)}")
    print("The branch points are the roots of this polynomial and the point at infinity.")
    print("\nOne root is clearly x=0. To find the others, we analyze the polynomial P(x)/x.")

    g_coeffs = {4: 8, 3: 4, 2: 4, 1: 1, 0: 8}
    print(f"Let g(x) = P(x)/x = {print_polynomial(g_coeffs, var='x')}")

    print("\nStep 2: Find the 2-adic valuations of the roots of g(x) using its Newton Polygon.")
    print("The Newton Polygon is constructed from the points (i, v2(a_i)) where a_i is the coefficient of x^i.")

    valuations = {i: v2(c) for i, c in g_coeffs.items()}
    print("Coefficients of g(x) and their 2-adic valuations:")
    for i in sorted(valuations.keys()):
        print(f"  Coefficient of x^{i} (a_{i}): {g_coeffs[i]:<2}, v2(a_{i}) = {valuations[i]}")

    print("\nAnalyzing the lower convex hull of the points (i, v2(a_i)):")
    print(" - Segment 1 connects (0, 3) and (1, 0). Slope = -3. Length = 1.")
    print("   This implies one root 'r1' with valuation v2(r1) = -(-3) = 3.")
    print(" - Segment 2 connects (1, 0) and (4, 3) (passing through (3,2)). Slope = 1. Length = 3.")
    print("   This implies three roots 'r2, r3, r4' with valuation v2(r_i) = -(1) = -1.")

    print("\nStep 3: Cluster the 6 branch points of the original curve.")
    print("The branch points are {0, r1, r2, r3, r4, infinity}.")
    print("We cluster them based on whether their 2-adic valuation is >= 0 (integral) or < 0 (non-integral).")

    # A point at infinity is treated as non-integral in this context.
    # Cluster 1: Integral roots (v2 >= 0)
    # The root 0 has v2 = inf. The root r1 has v2 = 3.
    size_B1 = 2
    # Cluster 2: Non-integral roots (v2 < 0)
    # The roots r2, r3, r4 have v2 = -1. The point at infinity also belongs to this cluster.
    size_B2 = 4
    
    print(f"\nCluster 1 (integral points): {{0, r1}}. Size = {size_B1}")
    print(f"Cluster 2 (non-integral points): {{r2, r3, r4, infinity}}. Size = {size_B2}")

    print("\nStep 4: Determine the stable model and count double points.")
    print("The theory of hyperelliptic curve reduction states that when the branch points partition into two sets of even size (here, 2 and 4),")
    print("the stable model consists of two components intersecting transversally at multiple points.")
    print("The structure is an elliptic curve (genus 1) and a rational curve (genus 0) intersecting at 2 points.")
    print("These intersection points are the double points (nodes) of the stable reduction.")

    num_double_points = 2
    
    # Final equation structure explanation
    # C_stable = C_elliptic U C_rational
    # |C_elliptic intersect C_rational| = 2
    # Genus check: g_total = g_elliptic + g_rational + num_intersections - 1
    # 2 = 1 + 0 + 2 - 1 = 2. This is correct.
    
    print("\nFinal calculation:")
    print("The stable model has 2 components which intersect at a number of points.")
    print("Number of intersections = Number of double points.")
    print(f"In this case, the number of double points is {num_double_points}.")

main()