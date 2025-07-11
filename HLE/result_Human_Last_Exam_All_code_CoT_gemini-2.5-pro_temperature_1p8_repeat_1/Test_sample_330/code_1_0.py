import numpy as np

def solve_and_print_square_vertices():
    """
    This function restores a square from four points given on its sides.
    It performs the calculation and prints the steps and final result.
    """
    # The coordinates of the four marked points
    p1 = np.array([0.3511, 0.2027])
    p2 = np.array([0.6753, 0.8303])
    p3 = np.array([-0.2845, 0.9905])
    p4 = np.array([-0.128, 0.2218])

    print("Step 1: Set up the geometric conditions.")
    print("Assuming the points are on consecutive sides, P1 and P3 are opposite, and P2 and P4 are opposite.")
    
    # Calculate vectors between opposite points
    # Vector P1->P3
    v_13 = p3 - p1
    # Vector P2->P4
    v_24 = p4 - p2

    # The condition that the side lengths are equal leads to a quadratic equation
    # for t = tan(theta), where theta is the orientation angle of the square.
    # The equation has the form: A*t^2 + B*t + C = 0
    a, b = v_13[0], v_13[1]
    d, e = v_24[0], v_24[1]

    A_coeff = a**2 - e**2
    B_coeff = -2 * (a*b + d*e)
    C_coeff = b**2 - d**2
    
    print("\nStep 2: Solve the central quadratic equation for t = tan(theta).")
    print("The final equation used to determine the square's orientation is:")
    print(f"({A_coeff:.4f}) * t^2 + ({B_coeff:.4f}) * t + ({C_coeff:.4f}) = 0")
    
    # Solve for t = tan(theta)
    delta = B_coeff**2 - 4*A_coeff*C_coeff
    # There are two solutions for t, corresponding to two possible squares. We'll pick one.
    t_solution = (-B_coeff + np.sqrt(delta)) / (2*A_coeff)
    
    theta = np.arctan(t_solution)
    c, s = np.cos(theta), np.sin(theta)

    print(f"\nUsing one solution, t = {t_solution:.4f}, we find the angle of the square's sides.")

    print("\nStep 3: Determine the equations for the four sides of the square (Ax + By + C = 0).")
    # Line AB (side 1, through P1): s*x - c*y - s*p1[0] + c*p1[1] = 0
    L1 = [s, -c, -s*p1[0] + c*p1[1]]
    # Line BC (side 2, through P2): c*x + s*y - c*p2[0] - s*p2[1] = 0
    L2 = [c, s, -c*p2[0] - s*p2[1]]
    # Line CD (side 3, through P3): s*x - c*y - s*p3[0] + c*p3[1] = 0
    L3 = [s, -c, -s*p3[0] + c*p3[1]]
    # Line DA (side 4, through P4): c*x + s*y - c*p4[0] - s*p4[1] = 0
    L4 = [c, s, -c*p4[0] - s*p4[1]]
    print("Side 1: {0:.4f}x + {1:.4f}y + {2:.4f} = 0".format(L1[0], L1[1], L1[2]))
    print("Side 2: {0:.4f}x + {1:.4f}y + {2:.4f} = 0".format(L2[0], L2[1], L2[2]))
    print("Side 3: {0:.4f}x + {1:.4f}y + {2:.4f} = 0".format(L3[0], L3[1], L3[2]))
    print("Side 4: {0:.4f}x + {1:.4f}y + {2:.4f} = 0".format(L4[0], L4[1], L4[2]))
    
    def intersect(line1, line2):
        A1, B1, C1 = line1; A2, B2, C2 = line2
        det = A1*B2 - A2*B1
        x = (B1*C2 - B2*C1) / det
        y = (A2*C1 - A1*C2) / det
        return np.array([x, y])

    # Find vertices by intersecting the lines
    V_A = intersect(L4, L1)
    V_B = intersect(L1, L2)
    V_C = intersect(L2, L3)
    V_D = intersect(L3, L4)
    
    vertices = [V_A, V_B, V_C, V_D]
    
    # Sort vertices by their x-coordinate
    vertices.sort(key=lambda v: v[0])
    
    print("\nStep 4: The coordinates of the vertices of the square, sorted by x-coordinate:")
    for v in vertices:
        print(f"({v[0]:.2f}, {v[1]:.2f})")

# Execute the main function
solve_and_print_square_vertices()