import math

def get_valuation(n, p):
    """Computes the p-adic valuation of an integer n."""
    if n == 0:
        return float('inf')
    valuation = 0
    while n % p == 0:
        valuation += 1
        n //= p
    return valuation

def newton_polygon_slopes(coeffs, p):
    """
    Computes the slopes of the Newton polygon for a polynomial with given coefficients.
    Coeffs should be a list [c_0, c_1, ..., c_n].
    """
    points = []
    print("Polynomial coefficients and their 2-adic valuations:")
    for i, c in enumerate(coeffs):
        if c != 0:
            val = get_valuation(c, p)
            points.append((i, val))
            print(f"  - Coeff of x^{i} (c_{i}): {c}, v_2(c_{i}) = {val}")
    
    # Sort points by x-coordinate
    points.sort()
    
    print("\nPoints for the Newton Polygon (degree, valuation):")
    print(f"  {points}")

    # Find the lower convex hull (Newton polygon)
    hull = [points[0]]
    for i in range(1, len(points)):
        while len(hull) >= 2:
            p1 = hull[-2]
            p2 = hull[-1]
            p3 = points[i]
            # Check if p3 is to the left of the line p1-p2 (cross product)
            # (p2[0]-p1[0])*(p3[1]-p1[1]) - (p2[1]-p1[1])*(p3[0]-p1[0])
            if (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]) >= 0:
                hull.pop()
            else:
                break
        hull.append(points[i])
        
    print("\nPoints on the lower convex hull (Newton Polygon):")
    print(f"  {hull}")

    slopes = []
    print("\nSegments of the Newton Polygon and their slopes:")
    for i in range(len(hull) - 1):
        p1 = hull[i]
        p2 = hull[i+1]
        slope = (p2[1] - p1[1]) / (p2[0] - p1[0])
        length = p2[0] - p1[0]
        print(f"  - Segment from {p1} to {p2}:")
        print(f"    - Slope = ({p2[1]} - {p1[1]}) / ({p2[0]} - {p1[0]}) = {slope}")
        print(f"    - Horizontal length = {length}")
        slopes.append({'slope': slope, 'length': length})
        
    return slopes

# Coefficients of P(x) = 8x^4 + 4x^3 + 4x^2 + x + 8
# The list is [c_0, c_1, c_2, c_3, c_4]
coeffs = [8, 1, 4, 4, 8]
p = 2

slopes = newton_polygon_slopes(coeffs, p)

print("\nValuations of the roots of P(x):")
for s in slopes:
    num_roots = s['length']
    valuation = -s['slope']
    print(f"  - There are {num_roots} root(s) with 2-adic valuation {valuation}.")

print("\n--- Analysis of the Stable Reduction ---")
print("The 6 branch points of the curve y^2=f(x) are the 5 roots of f(x) and the point at infinity.")
print("The valuations of the branch points are:")
print("  - Root x=0: v_2(0) = infinity")
print("  - Root from P(x) (alpha_1): v_2(alpha_1) = 3.0")
print("  - Three roots from P(x) (beta_i): v_2(beta_i) = -1.0")
print("  - Point at infinity: Effectively has valuation -infinity")

print("\nThese branch points form two clusters:")
print("  - Cluster S1 = {0, alpha_1}: The two points 'close to 0'. Size = 2.")
print("  - Cluster S2 = {beta_1, beta_2, beta_3, infinity}: The four points 'far from 0'. Size = 4.")

print("\nThe stable reduction has two components, C1 and C2, corresponding to these clusters.")
print("C1 is a double cover of the projective line branched at 2 points. Its genus is (2/2) - 1 = 0.")
print("C2 is a double cover of the projective line branched at 4 points. Its genus is (4/2) - 1 = 1.")

print("\nThe two components C1 (genus 0) and C2 (genus 1) must intersect to form a connected curve.")
print("Since the cover is 2-to-1 and unramified at the intersection point of the base components, C1 and C2 intersect at 2 points.")
print("These two intersection points are the ordinary double points (nodes) of the stable reduction.")

print("\nFinal calculation:")
num_double_points = 2
print(f"The number of double points is {num_double_points}.")
