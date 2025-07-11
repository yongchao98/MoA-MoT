import numpy as np

def verify_b():
    """
    This function verifies the construction used to answer part (b).
    """
    print("--- Verification for Part (b) ---")
    
    # Define vectors u and y in R^14. We only need the first 5 components.
    u = np.zeros(14)
    y = np.zeros(14)
    u[0:3] = [2, 2, 1]
    y[0:5] = [2, 2, 1, 6, 3]

    print(f"Vector u = {u}")
    print(f"Vector y = {y}")

    # Check properties of u
    u_dot_u = np.dot(u, u)
    print(f"\nChecking vector u:")
    print(f"u . u = {u_dot_u}")
    if u_dot_u == 9:
        print("PASS: u.u is 9 as required.")
    else:
        print(f"FAIL: u.u is not 9.")

    # Check properties of y
    y_dot_y = np.dot(y, y)
    print(f"\nChecking vector y:")
    print(f"y . y = {y_dot_y}")
    if y_dot_y == 54:
        print("PASS: y.y is 54, which is a multiple of 54.")
    else:
        print(f"FAIL: y.y is not 54.")

    # Define vector x = y/3
    x = y / 3.0
    print(f"\nVector x = y/3 = {x}")
    
    # Check properties of x
    x_dot_x = np.dot(x, x)
    print(f"\nChecking vector x:")
    print(f"x . x = {x_dot_x:.2f}")
    if abs(x_dot_x - 6.0) < 1e-9:
        print("PASS: x.x is 6, which is a multiple of 6.")
    else:
        print("FAIL: x.x is not 6.")

    # Check conditions for x to be in the lattice L defined by u
    y_dot_u = np.dot(y, u)
    print(f"\nChecking conditions for x = y/3 to be in L:")
    print(f"y . u = {y_dot_u}")
    if y_dot_u % 9 == 0:
        print("PASS: y.u is a multiple of 9.")
    else:
        print("FAIL: y.u is not a multiple of 9.")

    y_mod_3 = y % 3
    u_mod_3 = u % 3
    print(f"y mod 3 = {y_mod_3}")
    print(f"u mod 3 = {u_mod_3}")
    if np.array_equal(y_mod_3, u_mod_3):
        print("PASS: Congruence condition y = j*u (mod 3) is met (for j=1).")
    else:
        print("FAIL: Congruence condition is not met.")
    
    # Check for 3-primitivity
    print(f"\nChecking if x is 3-primitive (i.e., x/3 is not in L):")
    # A necessary condition for a vector w to be in L is that 3*w is an integer vector.
    # Let w = x/3. Then 3*w = x.
    is_integer_vector = np.all(np.equal(np.mod(x, 1), 0))
    print(f"Is x an integer vector? {is_integer_vector}")
    if not is_integer_vector:
        print("PASS: x is not an integer vector, so 3*(x/3) is not an integer vector.")
        print("This means x/3 cannot be in the lattice L, so x is 3-primitive.")
    else:
        print("FAIL: Could not prove 3-primitivity with this argument.")

verify_b()