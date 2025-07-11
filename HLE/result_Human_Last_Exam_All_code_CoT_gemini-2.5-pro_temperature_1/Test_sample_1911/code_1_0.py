def analyze_set_L_minus_point(z):
    """
    Analyzes the set L \ {z} and explains whether it can be a smooth manifold.
    L is the set {(x, y) | y = |x|}.
    z is a point (x, y).
    """
    x, y = z
    
    # Check if the point z is on the set L
    if y != abs(x):
        print(f"The point z = ({x}, {y}) is not in L, since y != |x|.")
        return

    print(f"The point z = ({x}, {y}) is in L.")

    # Case 1: The point z is the origin
    if x == 0 and y == 0:
        print("This point is the origin (0, 0).")
        print("The set L \\ {z} is L \\ {(0, 0)}.")
        print("This set consists of two disjoint open rays: {(t, t) | t > 0} and {(-t, t) | t > 0}.")
        print("Each ray is a smooth 1-dimensional manifold (diffeomorphic to the real line).")
        print("The disjoint union of smooth manifolds is a smooth manifold.")
        print("Therefore, L \\ {(0, 0)} can be given the structure of a smooth manifold.")
    
    # Case 2: The point z is not the origin
    else:
        print(f"This point is z = ({x}, {y}), which is not the origin.")
        print(f"The set L \\ {{z}} is L \\ {{({x}, {y})}}.")
        print("This set still contains the origin (0, 0), which forms a 'corner'.")
        print("A neighborhood of the origin in this set is not locally Euclidean.")
        print("Therefore, L \\ {z} cannot be given the structure of a smooth manifold.")

def main():
    """
    Demonstrates the logic for statement E.
    E: There exists a unique z in L such that L \ {z} can be given the structure of a smooth manifold.
    """
    print("--- Analysis for Statement E ---")
    
    # Test the unique point z = (0,0)
    z1_x, z1_y = 0, 0
    print(f"\nTesting the point z = ({z1_x}, {z1_y})...")
    analyze_set_L_minus_point((z1_x, z1_y))
    
    # Test an example of another point z != (0,0)
    z2_x, z2_y = 3, 3
    print(f"\nTesting the point z = ({z2_x}, {z2_y})...")
    analyze_set_L_minus_point((z2_x, z2_y))
    
    print("\n--- Conclusion for Statement E ---")
    print("The analysis shows that only when z = (0, 0) is removed can the remaining set be a smooth manifold.")
    print("This means the point z = (0, 0) is unique.")
    print("Therefore, statement E is TRUE.")

if __name__ == "__main__":
    main()
