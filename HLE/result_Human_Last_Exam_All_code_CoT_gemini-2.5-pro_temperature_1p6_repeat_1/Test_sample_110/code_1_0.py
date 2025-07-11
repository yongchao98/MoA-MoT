import math

def solve_robot_arm_distance():
    """
    Calculates the minimum distance between a robot arm's shoulder and finger.
    """
    # Segment lengths
    L1 = 40  # Shoulder to Elbow
    L2 = 28  # Elbow to Wrist
    L3 = 15  # Wrist to Hand
    L4 = 10  # Hand to Finger

    print(f"Segment lengths: L1={L1}, L2={L2}, L3={L3}, L4={L4}\n")

    # --- Step 1: Determine joint angles from constraints ---
    # The primary constraint is the 1cm minimum distance between non-adjacent segments.
    # We model this as: L_separator * sin(d_joint) = 1 cm
    # This gives us the minimum angle for maximum folding.

    # Angle d_E (Elbow) is constrained by L1-L3 distance, separated by L2.
    sin_d_E = 1 / L2
    d_E = math.asin(sin_d_E)
    cos_d_E = math.cos(d_E)
    print(f"Constraint dist(L1, L3) >= 1cm, separated by L2={L2}cm.")
    print(f"  - sin(d_E) = 1 / {L2} = {sin_d_E:.4f}")
    print(f"  - d_E = {math.degrees(d_E):.4f} degrees\n")

    # Angle d_W (Wrist) is constrained by L2-L4 distance, separated by L3.
    sin_d_W = 1 / L3
    d_W = math.asin(sin_d_W)
    cos_d_W = math.cos(d_W)
    print(f"Constraint dist(L2, L4) >= 1cm, separated by L3={L3}cm.")
    print(f"  - sin(d_W) = 1 / {L3} = {sin_d_W:.4f}")
    print(f"  - d_W = {math.degrees(d_W):.4f} degrees\n")

    # Angle d_H (Hand) is unconstrained by a subsequent segment, so we assume d_H=0 for maximum folding.
    d_H = 0
    print("Assuming d_H = 0 for maximum folding.\n")

    # --- Step 2: Calculate coordinates of the fingertip (F) ---
    # We place the shoulder S at (0,0) and the first segment L1 along the x-axis.
    # The arm folds back in a zig-zag pattern.
    # The final coordinates (Fx, Fy) are the sum of the segment vectors.
    # Fx = L1 - L2*cos(d_E) + L3*cos(d_E+d_W) - L4*cos(d_E+d_W+d_H)
    # Fy = L2*sin(d_E) - L3*sin(d_E+d_W) + L4*sin(d_E+d_W+d_H)

    # Since d_H=0, the equations simplify:
    # Fx = L1 - L2*cos(d_E) + (L3-L4)*cos(d_E+d_W)
    # Fy = L2*sin(d_E) - (L3-L4)*sin(d_E+d_W)
    
    # Pre-calculate trigonometric values for the combined angle d_E+d_W
    s_EW = math.sin(d_E + d_W)
    c_EW = math.cos(d_E + d_W)
    
    term_L2_x = L2 * cos_d_E
    term_L2_y = L2 * sin_d_E # This is 28 * (1/28) = 1

    term_L3_minus_L4_x = (L3 - L4) * c_EW
    term_L3_minus_L4_y = (L3 - L4) * s_EW

    Fx = L1 - term_L2_x + term_L3_minus_L4_x
    Fy = term_L2_y - term_L3_minus_L4_y
    
    print("Calculating final coordinates of the finger (F):")
    print(f"  Fx = L1 - L2*cos(d_E) + (L3-L4)*cos(d_E+d_W)")
    print(f"     = {L1} - {L2}*{cos_d_E:.4f} + ({L3}-{L4})*{c_EW:.4f}")
    print(f"     = {L1} - {term_L2_x:.4f} + {term_L3_minus_L4_x:.4f} = {Fx:.4f}\n")

    print(f"  Fy = L2*sin(d_E) - (L3-L4)*sin(d_E+d_W)")
    print(f"     = {L2}*{sin_d_E:.4f} - ({L3}-{L4})*{s_EW:.4f}")
    print(f"     = {term_L2_y:.4f} - {term_L3_minus_L4_y:.4f} = {Fy:.4f}\n")

    # --- Step 3: Calculate the final distance ---
    distance = math.sqrt(Fx**2 + Fy**2)
    print("Calculating final distance from Shoulder S(0,0) to Finger F(Fx,Fy):")
    print(f"  Distance = sqrt(Fx^2 + Fy^2)")
    print(f"           = sqrt({Fx:.4f}^2 + {Fy:.4f}^2)")
    print(f"           = {distance:.4f} cm")

solve_robot_arm_distance()