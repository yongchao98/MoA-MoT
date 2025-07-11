import numpy as np

def explain_l1_geometry():
    """
    This script explains and calculates the angle between the solar wind flow
    and the interplanetary magnetic field (IMF) at the L1 point.
    """

    print("--- Background ---")
    print("Your question: Why is normalized magnetic helicity calculated using components")
    print("perpendicular to the radial direction when AIC waves propagate along the local magnetic field?")
    print("\nThis script demonstrates the core reason: The local magnetic field is not radial.")
    print("At L1, the magnetic field follows the Parker Spiral, making an angle with the radial flow.")

    # 1. Define the vectors.
    # The solar wind flows radially from the Sun. We can represent this with a unit vector.
    radial_direction = np.array([1.0, 0.0, 0.0])

    # The Parker Spiral angle at 1 AU is typically ~45 degrees from the radial direction
    # in the ecliptic plane.
    parker_spiral_angle_deg = 45.0
    parker_spiral_angle_rad = np.deg2rad(parker_spiral_angle_deg)
    
    # We define the magnetic field direction vector in the same plane for simplicity.
    b_field_direction = np.array([np.cos(parker_spiral_angle_rad), np.sin(parker_spiral_angle_rad), 0.0])

    # 2. Calculate the angle between the two vectors.
    # Formula: angle = arccos( (A · B) / (|A| * |B|) )
    dot_product = np.dot(radial_direction, b_field_direction)
    mag_radial = np.linalg.norm(radial_direction)
    mag_b_field = np.linalg.norm(b_field_direction)
    
    # We clip the argument to arccos to avoid floating point inaccuracies if it's slightly > 1.0
    cos_theta = np.clip(dot_product / (mag_radial * mag_b_field), -1.0, 1.0)
    angle_rad_result = np.arccos(cos_theta)
    angle_deg_result = np.rad2deg(angle_rad_result)

    # 3. Print the step-by-step calculation.
    print("\n--- Step-by-Step Calculation ---")
    print("We calculate the angle between the radial flow and the magnetic field.")
    print("\nEquation: angle = arccos( (Radial_Dir · B_Field_Dir) / (|Radial_Dir| * |B_Field_Dir|) )")

    print(f"\n1. Radial Direction Vector = [{radial_direction[0]}, {radial_direction[1]}, {radial_direction[2]}]")
    print(f"2. B-Field Direction Vector = [{b_field_direction[0]:.4f}, {b_field_direction[1]:.4f}, {b_field_direction[2]:.4f}]")
    
    print(f"\n3. Dot Product (Radial_Dir · B_Field_Dir):")
    print(f"   ({radial_direction[0]} * {b_field_direction[0]:.4f}) + ({radial_direction[1]} * {b_field_direction[1]:.4f}) + ({radial_direction[2]} * {b_field_direction[2]:.4f}) = {dot_product:.4f}")

    print(f"\n4. Magnitude of Radial_Dir: |Radial_Dir| = {mag_radial}")
    print(f"5. Magnitude of B_Field_Dir: |B_Field_Dir| = {mag_b_field:.4f}")
    
    print(f"\n6. Putting it all together:")
    print(f"   angle = arccos( {dot_product:.4f} / ({mag_radial} * {mag_b_field:.4f}) )")
    print(f"   angle = arccos( {cos_theta:.4f} )")
    print(f"   angle = {angle_deg_result:.2f} degrees")

    print("\n--- Conclusion ---")
    print(f"The angle is {angle_deg_result:.2f} degrees, not zero. Therefore, the magnetic field is not radial.")
    print("Calculating helicity in the plane perpendicular to the radial direction is an approximation")
    print("justified by the Taylor 'frozen-in flow' Hypothesis, where the solar wind's radial flow")
    print("dominates the measurement.")

# Run the explanation
explain_l1_geometry()