import math

def get_valuation(n):
    """Computes the 2-adic valuation of an integer."""
    if n == 0:
        return float('inf')
    valuation = 0
    while n % 2 == 0:
        valuation += 1
        n //= 2
    return valuation

def main():
    """
    Calculates the thickness of the double point of the stable reduction of the curve
    z^2 = 2*x^5 + 2*x^3 + 1.
    """
    print("Step 1: Define the polynomial from the curve equation.")
    # Original polynomial P(x) = 2*x^5 + 2*x^3 + 1
    print("The curve is z^2 = 2*x^5 + 2*x^3 + 1.")
    print("\nStep 2: Transform the polynomial to be monic.")
    print("We use the substitution x = y/2.")
    print("The equation becomes 16*z^2 = y^5 + 4*y^3 + 16.")
    print("Let Z = 4z, the new model is Z^2 = f(y), where f(y) = y^5 + 4*y^3 + 16.")
    
    # Coefficients of f(y)
    coeffs = {0: 16, 1: 0, 2: 0, 3: 4, 4: 0, 5: 1}
    
    print("\nStep 3: Determine the Newton Polygon and root valuations.")
    points = []
    for i, coeff in coeffs.items():
        if coeff != 0:
            points.append((i, get_valuation(coeff)))
    
    # Sort points by degree
    points.sort()
    print(f"The points for the Newton Polygon are (degree, valuation): {points}")
    
    # The points for the lower convex hull are (0,4), (3,2), (5,0)
    hull_points = [(0, 4), (3, 2), (5, 0)]
    print(f"The vertices of the lower convex hull are: {hull_points}")
    
    # Calculate slopes and lengths
    slopes = []
    for i in range(len(hull_points) - 1):
        p1 = hull_points[i]
        p2 = hull_points[i+1]
        slope = (p2[1] - p1[1]) / (p2[0] - p1[0])
        length = p2[0] - p1[0]
        slopes.append({'slope': slope, 'length': length})
        
    print("The slopes of the Newton Polygon segments are:")
    for s in slopes:
        print(f"  Slope: {s['slope']}, horizontal length: {s['length']}")
        
    root_valuations = {}
    print("The valuations of the roots are the negatives of the slopes:")
    for s in slopes:
        valuation = -s['slope']
        num_roots = s['length']
        root_valuations[valuation] = num_roots
        print(f"  Valuation: {valuation}, number of roots: {num_roots}")
        
    print("\nStep 4: Identify the clusters of roots.")
    # The valuations are 2/3 and 1. These form two distinct clusters.
    cluster1_val = 2/3
    cluster1_size = root_valuations[cluster1_val]
    cluster2_val = 1.0
    cluster2_size = root_valuations[cluster2_val]
    
    print(f"Cluster C1: {cluster1_size} roots with valuation {cluster1_val:.2f}")
    print(f"Cluster C2: {cluster2_size} roots with valuation {cluster2_val}")
    
    print("\nStep 5: Calculate the thickness of the double point.")
    # The potential pi(C1, C2) = |C1|*|C2|*min(v1, v2)
    min_val = min(cluster1_val, cluster2_val)
    potential = cluster1_size * cluster2_size * min_val
    
    print("The potential between the clusters is pi = |C1| * |C2| * min(v1, v2)")
    print(f"pi = {cluster1_size} * {cluster2_size} * {min_val:.2f} = {potential}")
    
    # The thickness delta = floor(potential / 2)
    thickness = math.floor(potential / 2)
    
    print("The thickness of the double point is given by delta = floor(pi / 2)")
    print(f"delta = floor({potential} / 2)")
    print(f"delta = floor({potential/2})")
    print(f"delta = {thickness}")
    
    print("\nFinal Answer:")
    print(f"The thickness of the double point is {thickness}.")

if __name__ == "__main__":
    main()
