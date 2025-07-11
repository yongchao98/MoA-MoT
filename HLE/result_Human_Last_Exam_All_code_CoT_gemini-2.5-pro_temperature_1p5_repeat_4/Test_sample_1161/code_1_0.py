import numpy as np

def check_and_print_coverage(guards, point):
    """
    Checks if a point is seen by any guard and prints the results.
    A point p is seen by a guard g if the dot product p·g >= 1.
    """
    print(f"Testing point P = {np.round(point, 4)}")
    print(f"Distance of P from origin: {np.linalg.norm(point):.4f} > 1")
    
    is_seen = False
    print("\n--- Checking Coverage ---")
    for i, guard in enumerate(guards):
        dot_product = np.dot(point, guard)
        if dot_product >= 1:
            is_seen = True
            print(f"P · G_{i+1} = {dot_product:.4f} >= 1  (Point is SEEN by this guard)")
        else:
            print(f"P · G_{i+1} = {dot_product:.4f} < 1")

    print("--- Result ---")
    if is_seen:
        print("The point P IS SEEN by at least one guard.")
    else:
        print("The point P IS NOT SEEN by any guard.")
    print("="*40 + "\n")

def run_demonstration():
    """
    Runs demonstrations for N=4 and N=6 guards.
    """
    # Case 1: N=4 guards at the vertices of an inscribed regular tetrahedron.
    print("CASE 1: N=4 Guards (Tetrahedral arrangement)")
    # Vertices of a regular tetrahedron on the unit sphere
    g1 = np.array([np.sqrt(8/9), 0, -1/3])
    g2 = np.array([-np.sqrt(2/9), np.sqrt(2/3), -1/3])
    g3 = np.array([-np.sqrt(2/9), -np.sqrt(2/3), -1/3])
    g4 = np.array([0, 0, 1])
    guards_4 = [g1, g2, g3, g4]

    # An unseen "hole" exists in the direction of the normal to a face.
    # The normal to the face (g1, g2, g3) is in the direction of -g4.
    # The unseen region extends up to r=3 in this direction.
    # Let's test a point at r=2 in this direction.
    r = 2.0
    hole_direction_4 = -g4
    test_point_4 = r * hole_direction_4
    check_and_print_coverage(guards_4, test_point_4)

    # Case 2: N=6 guards at the vertices of an inscribed regular octahedron.
    print("CASE 2: N=6 Guards (Octahedral arrangement)")
    # Vertices of a regular octahedron on the unit sphere
    guards_6 = [
        np.array([1, 0, 0]), np.array([-1, 0, 0]),
        np.array([0, 1, 0]), np.array([0, -1, 0]),
        np.array([0, 0, 1]), np.array([0, 0, -1])
    ]

    # An unseen "hole" exists in the direction of a face normal (1,1,1).
    # The unseen region extends up to r=sqrt(3) (~1.732) in this direction.
    # Let's test a point at r=1.5 in this direction.
    r = 1.5
    hole_direction_6 = np.array([1, 1, 1]) / np.sqrt(3)
    test_point_6 = r * hole_direction_6
    check_and_print_coverage(guards_6, test_point_6)
    
if __name__ == '__main__':
    run_demonstration()
