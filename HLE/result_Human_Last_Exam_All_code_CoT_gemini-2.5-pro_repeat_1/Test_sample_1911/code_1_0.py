import math

def demonstrate_unboundedness():
    """
    Demonstrates that the set L = {(x, |x|)} is not bounded, and thus not compact.
    
    A key topological property is compactness. For a set to be diffeomorphic
    to an n-sphere (S^n), it must be compact. In R^2, a necessary condition
    for compactness is being bounded. This function shows L is not bounded
    by generating points on L that are increasingly far from the origin.
    """
    print("Statement C claims L can be diffeomorphic to S^n for any n.")
    print("This is false because S^n is compact, but L is not.")
    print("A compact set in R^2 must be bounded. Let's show L is unbounded.")
    print("-" * 60)
    print("We will generate points p_k = (k, |k|) on L for increasing k.")
    print("Then we calculate their distance (norm) from the origin (0,0).")
    print("-" * 60)

    for k in [1, 10, 100, 1000, 10000]:
        # A point on L is (x, y) where y = |x|
        point_x = float(k)
        point_y = abs(point_x)
        
        # The distance from the origin is sqrt(x^2 + y^2)
        # distance = sqrt(k^2 + k^2) = sqrt(2 * k^2) = k * sqrt(2)
        distance = math.sqrt(point_x**2 + point_y**2)
        
        print(f"For k = {k}:")
        print(f"  Point on L: ({point_x:.0f}, {point_y:.0f})")
        print(f"  Distance from origin = sqrt({point_x:.0f}^2 + {point_y:.0f}^2) = {distance:.2f}")

    print("-" * 60)
    print("As k increases, the distance from the origin grows without bound.")
    print("Therefore, the set L is not bounded and cannot be compact.")
    print("Since L is not compact and S^n is, they cannot be homeomorphic,")
    print("let alone diffeomorphic. Thus, statement C is false.")

demonstrate_unboundedness()
