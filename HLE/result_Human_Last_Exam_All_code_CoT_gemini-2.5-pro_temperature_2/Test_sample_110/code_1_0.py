import math

def calculate_robot_arm_distance():
    """
    Calculates the minimum distance between the robot's finger and shoulder.
    """
    # Segment lengths in cm
    l1 = 40  # Shoulder to Elbow
    l2 = 28  # Elbow to Wrist
    l3 = 15  # Wrist to Hand
    l4 = 10  # Hand (finger)

    # Clearance constraint, interpreted as vertical separation between segments
    h = 3.5

    # Step-by-step calculation of the coordinates for each joint
    # P0 (Shoulder) is at the origin (0, 0)
    p0_x, p0_y = 0, 0

    # P1 (Elbow) is placed along the x-axis
    p1_x, p1_y = l1, 0
    print(f"Shoulder P0: ({p0_x}, {p0_y})")
    print(f"Segment L1 (length {l1})")
    print(f"Elbow P1: ({p1_x}, {p1_y})")
    print("-" * 20)

    # P2 (Wrist) position calculation
    # Folds back from P1. The segment forms a right triangle with sides (p1_x - p2_x), p2_y, and hypotenuse l2.
    # We set the vertical position of P2 based on the separation h.
    p2_y = h
    # Calculate the horizontal position p2_x using Pythagorean theorem
    delta_x2 = math.sqrt(l2**2 - (p2_y - p1_y)**2)
    p2_x = p1_x - delta_x2  # Subtract because it folds back
    print(f"Segment L2 (length {l2})")
    print(f"Wrist P2: (x={p2_x:.2f}, y={p2_y})")
    print(f"Equation for P2: ({p1_x} - x)^2 + ({p2_y} - {p1_y})^2 = {l2}^2")
    print(f"({p1_x} - {p2_x:.2f})^2 + ({p2_y} - {p1_y})^2 = {delta_x2:.2f}^2 + {p2_y**2} = {delta_x2**2+p2_y**2:.2f} (approx {l2**2})")
    print("-" * 20)
    
    # P3 (Hand) position calculation
    # Folds back from P2. Vertical position is another 'h' step away.
    p3_y = p2_y + h
    # Calculate horizontal position p3_x
    delta_x3 = math.sqrt(l3**2 - (p3_y - p2_y)**2)
    p3_x = p2_x - delta_x3
    print(f"Segment L3 (length {l3})")
    print(f"Hand Base P3: (x={p3_x:.2f}, y={p3_y})")
    print(f"Equation for P3: ({p2_x:.2f} - x)^2 + ({p3_y} - {p2_y})^2 = {l3}^2")
    print(f"({p2_x:.2f} - {p3_x:.2f})^2 + ({p3_y} - {p2_y})^2 = {delta_x3:.2f}^2 + {h**2} = {delta_x3**2+h**2:.2f} (approx {l3**2})")
    print("-" * 20)

    # P4 (Finger Tip) position calculation
    # Folds forward from P3.
    p4_y = p3_y + h
    # Calculate horizontal position p4_x
    delta_x4 = math.sqrt(l4**2 - (p4_y - p3_y)**2)
    p4_x = p3_x + delta_x4  # Add because it folds forward
    print(f"Segment L4 (length {l4})")
    print(f"Finger Tip P4: (x={p4_x:.2f}, y={p4_y})")
    print(f"Equation for P4: (x - {p3_x:.2f})^2 + ({p4_y} - {p3_y})^2 = {l4}^2")
    print(f"({p4_x:.2f} - {p3_x:.2f})^2 + ({p4_y} - {p3_y})^2 = {delta_x4:.2f}^2 + {h**2} = {delta_x4**2+h**2:.2f} (approx {l4**2})")
    print("-" * 20)

    # Final distance calculation
    final_distance = math.sqrt(p4_x**2 + p4_y**2)

    # Print the final result
    print("Final Calculation:")
    print(f"Distance = sqrt(P4_x^2 + P4_y^2)")
    print(f"Distance = sqrt({p4_x:.2f}^2 + {p4_y}^2)")
    print(f"Distance = sqrt({p4_x**2:.2f} + {p4_y**2:.2f})")
    print(f"Distance = sqrt({p4_x**2 + p4_y**2:.2f})")
    print(f"\nFinal calculated distance: {final_distance:.2f} cm")


calculate_robot_arm_distance()