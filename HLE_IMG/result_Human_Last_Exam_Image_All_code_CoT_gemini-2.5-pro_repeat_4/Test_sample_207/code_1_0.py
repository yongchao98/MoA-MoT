import sympy

def solve_robotics_inertia():
    """
    This function defines the symbolic variables and prints the expression for
    the M_32 element of the RPR robot's inertia matrix.
    """
    # Define the symbolic variables
    m3 = sympy.Symbol('m_3')       # Mass of link 3
    d_c3 = sympy.Symbol('d_c3')     # CoM distance for link 3
    q3 = sympy.Symbol('q3')         # Joint variable for joint 3

    # Construct the expression for M_32
    M_32 = m3 * d_c3 * sympy.sin(q3)

    # Print the result in a formatted equation
    print(f"The expression for the entry M_32 of the robot inertia matrix is:")
    print(f"M_32 = {M_32}")

solve_robotics_inertia()
