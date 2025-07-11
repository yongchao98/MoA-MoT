import numpy as np

def generate_cantor_set_points(iterations):
    """
    Generates a finite set of points approximating the Cantor set.
    The Cantor set is constructed by repeatedly removing the middle third.
    The endpoints of the remaining intervals are points in the Cantor set.
    """
    endpoints = {0.0, 1.0}
    intervals = [(0.0, 1.0)]
    
    for _ in range(iterations):
        new_intervals = []
        for start, end in intervals:
            length = end - start
            p1 = start + length / 3.0
            p2 = start + 2 * length / 3.0
            endpoints.add(p1)
            endpoints.add(p2)
            new_intervals.append((start, p1))
            new_intervals.append((p2, end))
        intervals = new_intervals
        
    return sorted(list(endpoints))

def demonstrate_non_open_component():
    """
    Demonstrates that components in C x (0,1) are not open.
    """
    print("This demonstration uses an analogy: the space X = C x [0,1], where C is the Cantor set.")
    print("The open set is V = C x (0,1).")
    print("The components of V are the vertical segments {c} x (0,1) for each c in C.")
    print("We will show that these components are not open.\n")

    # 1. Generate an approximation of the Cantor set
    # More iterations give a better approximation of the Cantor set's properties
    iterations = 4
    cantor_points = generate_cantor_set_points(iterations)
    
    # 2. Select a point 'c' in the Cantor set to define a component Kc
    # We pick a point that is not 0 or 1 to better illustrate.
    c = cantor_points[1] 
    print(f"Let's focus on the component K_c for c = {c:.4f}.")
    print(f"This component is the set of points {{ (x, y) | x = {c:.4f}, 0 < y < 1 }}.\n")

    # 3. Pick a point 'p' in this component
    p_x, p_y = c, 0.5
    print(f"Consider a point p = ({p_x:.4f}, {p_y:.2f}) within this component K_c.\n")

    # 4. Define a small neighborhood around p
    # In the space X = C x [0,1], a neighborhood is defined by an epsilon in both directions.
    epsilon = 0.02
    print(f"Let's define an open neighborhood N around p with radius epsilon = {epsilon:.4f}.")
    print(f"This neighborhood N contains all points (x, y) in X such that:")
    print(f"  {p_x - epsilon:.4f} < x < {p_x + epsilon:.4f}")
    print(f"  {p_y - epsilon:.2f} < y < {p_y + epsilon:.2f}\n")

    # 5. Show that this neighborhood contains points from other components
    # Find a point c' in the Cantor set that is close to c
    c_prime = -1
    for point in cantor_points:
        if 0 < abs(point - c) < epsilon:
            c_prime = point
            break
            
    if c_prime != -1:
        print(f"The Cantor set is 'dust-like'. Near our chosen c = {c:.4f}, there are other Cantor set points.")
        print(f"For instance, the point c' = {c_prime:.4f} is also in the Cantor set and is inside the x-range of our neighborhood N.")
        
        # This point p_prime is in the neighborhood N
        p_prime_x, p_prime_y = c_prime, 0.5
        
        print(f"\nNow consider the point p' = ({p_prime_x:.4f}, {p_prime_y:.2f}).")
        print(f"This point p' lies within the neighborhood N of p.")
        print(f"However, p' belongs to a different component, K_c' (where c' = {c_prime:.4f}).")
        print("It does not belong to our original component K_c.\n")
        
        print("Conclusion: The neighborhood N of p is not fully contained within the component K_c.")
        print("Since this is true for any point p and any neighborhood N, the component K_c is not an open set.\n")
    else:
        print("Could not find a nearby point with this epsilon. Try more iterations.")

    print("This logic applies to every component corresponding to every point in the Cantor set.")
    print("The number of points in the Cantor set is c (the cardinality of the continuum).")
    print("Therefore, the open set V = C x (0,1) has c non-open components.")
    print("\nThe largest possible number of non-open components is c.")

demonstrate_non_open_component()