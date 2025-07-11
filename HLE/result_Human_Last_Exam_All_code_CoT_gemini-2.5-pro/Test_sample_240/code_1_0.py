import math

def demonstrate_coupling_vs_angle():
    """
    This script demonstrates how the inductive coupling between a probe
    and a coil is dependent on their relative orientation.
    """
    print("This script demonstrates how geometric orientation affects coil coupling.")
    print("A coupling factor of 1.0 is maximum, and 0.0 is no coupling (a 'geometric null').\n")

    # A list of angles to test, in degrees.
    # 0 degrees: Probe field is parallel to the coil's sensitive axis (max coupling).
    # 90 degrees: Probe field is perpendicular (orthogonal) to the coil's axis (zero coupling).
    test_angles = [0, 30, 45, 60, 90]

    for angle_deg in test_angles:
        # Convert angle from degrees to radians for the math.cos function.
        angle_rad = math.radians(angle_deg)

        # The coupling strength is proportional to the absolute value of the cosine of the angle.
        coupling_factor = abs(math.cos(angle_rad))

        # Output the "equation" and the result for clarity.
        print(f"For an angle of {angle_deg} degrees:")
        # The vertical bars | | denote the absolute value.
        print(f"  Coupling Equation: Factor = |cos({angle_deg} degrees)|")
        print(f"  Resulting Coupling Factor = {coupling_factor:.4f}\n")

    print("As the results show, when the orientation reaches 90 degrees (orthogonal),")
    print("the coupling factor becomes 0. In this state, the coil's resonance would be invisible to the probe.")

if __name__ == '__main__':
    demonstrate_coupling_vs_angle()