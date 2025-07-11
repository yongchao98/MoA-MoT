import numpy as np

def check_concyclic():
    """
    Checks if the 5 leg points of the chair are concyclic.
    """
    p1 = np.array([0., 0.])
    p2 = np.array([2., 0.])
    p3 = np.array([2., 2.])
    p4 = np.array([0., 2.])
    p5 = np.array([1., 4.])

    points = {'P1': p1, 'P2': p2, 'P3': p3, 'P4': p4, 'P5': p5}
    
    print("Step 1: Define the five leg positions.")
    for name, pos in points.items():
        print(f"  {name}: {pos}")
    print("-" * 30)

    # We will use three non-collinear points P1, P2, P4 to define a circle.
    # The general equation of a circle is x^2 + y^2 + Dx + Ey + F = 0.
    # We solve a system of linear equations for D, E, F.
    
    # Matrix A for the system A * [D, E, F]' = B
    A = np.array([
        [p1[0], p1[1], 1],
        [p2[0], p2[1], 1],
        [p4[0], p4[1], 1]
    ])
    
    # Vector B
    B = -np.array([
        p1[0]**2 + p1[1]**2,
        p2[0]**2 + p2[1]**2,
        p4[0]**2 + p4[1]**2
    ])
    
    try:
        # Solve for D, E, F
        coeffs = np.linalg.solve(A, B)
        D, E, F = coeffs
        
        # Deriving h, k, r^2 from D, E, F
        # (x-h)^2 + (y-k)^2 = r^2  is x^2 - 2hx + y^2 - 2ky + h^2+k^2-r^2 = 0
        h = -D / 2
        k = -E / 2
        r_squared = h**2 + k**2 - F
        
        print("Step 2: Find the equation of the circle passing through P1, P2, and P4.")
        print(f"The equation is (x - {h})^2 + (y - {k})^2 = {r_squared}")
        print("-" * 30)

    except np.linalg.LinAlgError:
        print("The points P1, P2, and P4 are collinear and do not define a unique circle.")
        return

    print("Step 3: Test if the other points (P3 and P5) lie on this circle.")
    # Equation of the circle: (x-h)^2 + (y-k)^2 - r^2 = 0
    # A point is on the circle if this expression evaluates to 0.

    # Test P3
    val_p3 = (p3[0] - h)**2 + (p3[1] - k)**2 - r_squared
    print(f"Testing P3(2,2): (2 - {h})^2 + (2 - {k})^2 - {r_squared} = {val_p3:.2f}")
    if np.isclose(val_p3, 0):
        print("  > Result: P3 lies on the circle.")
    else:
        print("  > Result: P3 does NOT lie on the circle.")

    # Test P5
    val_p5 = (p5[0] - h)**2 + (p5[1] - k)**2 - r_squared
    print(f"Testing P5(1,4): (1 - {h})^2 + (4 - {k})^2 - {r_squared} = {val_p5:.2f}")
    if np.isclose(val_p5, 0):
        print("  > Result: P5 lies on the circle.")
    else:
        print("  > Result: P5 does NOT lie on the circle.")
    print("-" * 30)
    
    print("Conclusion: The five points are not concyclic.")
    print("This implies the chair cannot be placed on a perfect sphere, and for a generic uneven surface, the placement problem is over-constrained, meaning a solution is not guaranteed.")


if __name__ == '__main__':
    check_concyclic()
