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

def find_newton_polygon_slopes(coeffs, p):
    """
    Finds the slopes of the Newton Polygon for a polynomial.
    coeffs: A list of coefficients [c_0, c_1, ..., c_n].
    """
    points = []
    for i, c in enumerate(coeffs):
        v = get_valuation(c, p)
        if v != float('inf'):
            points.append((i, v))
    
    # Sort points by x-coordinate to be safe
    points.sort()

    # Find lower convex hull (Newton Polygon)
    if not points:
        return []

    hull = [points[0]]
    for i in range(1, len(points)):
        while len(hull) >= 2:
            p1 = hull[-2]
            p2 = hull[-1]
            p3 = points[i]
            # Check if p2 is "above" the line p1-p3 (cross product)
            # (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x)
            # We want points on or below the line, so we keep them if the cross product is >= 0
            if (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]) >= 0:
                hull.pop()
            else:
                break
        hull.append(points[i])
        
    slopes = []
    for i in range(len(hull) - 1):
        p1 = hull[i]
        p2 = hull[i+1]
        slope = (p2[1] - p1[1]) / (p2[0] - p1[0])
        width = p2[0] - p1[0]
        slopes.append({'slope': slope, 'width': width})
    
    return slopes

def main():
    """
    Main function to solve the problem.
    """
    print("Analyzing the curve y^2 = 8*x^5 + 4*x^4 + 4*x^3 + 1*x^2 + 8*x")
    print("This is a genus 2 hyperelliptic curve.")
    print("The number of double points in the stable reduction depends on the 2-adic valuations of the branch points.")
    print("The branch points are the roots of f(x) = 8*x^5 + 4*x^4 + 4*x^3 + x^2 + 8*x and the point at infinity.")
    print("\nFactoring f(x) = x * g(x), where g(x) = 8*x^4 + 4*x^3 + 4*x^2 + x + 8.")
    print("One branch point is x=0. The others are the four roots of g(x) and infinity.")
    
    # Coefficients of g(x) = 8x^4 + 4x^3 + 4x^2 + 1x + 8
    coeffs_g = [8, 1, 4, 4, 8] # [c0, c1, c2, c3, c4]
    p = 2
    
    slopes_info = find_newton_polygon_slopes(coeffs_g, p)
    
    print("\nTo find the valuations of the roots of g(x), we analyze its 2-adic Newton Polygon.")
    print(f"The coefficients of g(x) are {coeffs_g[::-1]} for powers x^4 down to x^0.")
    
    print("\nThe Newton Polygon has segments with the following slopes and widths:")
    for s in slopes_info:
        print(f"  - Slope: {s['slope']}, Width: {s['width']}")
    
    print("\nThe valuations of the roots are the negative of the slopes.")
    print("Therefore, the valuations of the roots of g(x) are:")
    for s in slopes_info:
        for _ in range(s['width']):
            print(f"  - Valuation: {-s['slope']}")
    
    print("\nSo we have one root with valuation 3 and three roots with valuation -1.")
    print("This suggests partitioning the 6 branch points {0, infinity, r1, r2, r3, r4} into two sets of three:")
    print("  S1 = {roots with valuation -1}")
    print("  S2 = {0, infinity, root with valuation 3}")
    print("For a genus 2 curve, this partition implies the stable reduction consists of two rational curves intersecting.")
    print("The number of intersections (double points), delta, is given by the formula g = g1 + g2 + delta - 1.")
    g = 2
    g1 = 0
    g2 = 0
    delta = g - g1 - g2 + 1
    print(f"With g = {g}, g1 = {g1}, g2 = {g2}, the formula is: {g} = {g1} + {g2} + delta - 1.")
    print(f"Solving for delta gives delta = {delta}.")

main()