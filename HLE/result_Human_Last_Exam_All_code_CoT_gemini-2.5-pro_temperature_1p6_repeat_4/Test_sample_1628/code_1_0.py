import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates and reports the shear center location for an asymmetric channel section.
    """
    # 1. Define Section Geometry (in mm)
    # An asymmetric channel with a taller web and wider top flange
    h = 150.0  # Height of the web
    tw = 6.0   # Thickness of the web
    b1 = 100.0 # Width of the top flange
    t1 = 10.0  # Thickness of the top flange
    b2 = 75.0  # Width of the bottom flange
    t2 = 8.0   # Thickness of the bottom flange

    print("--- Input Section Dimensions (mm) ---")
    print(f"Web (h x tw):      {h} x {tw}")
    print(f"Top Flange (b1 x t1): {b1} x {t1}")
    print(f"Bottom Flange (b2 x t2): {b2} x {t2}\n")

    # 2. Establish Coordinate System and Component Properties
    # Origin is at the top-left outer corner of the section.
    # y is positive downwards, x is positive to the right.

    # Top Flange
    A1 = b1 * t1
    y1 = t1 / 2
    x1 = b1 / 2

    # Web
    Aw = h * tw
    yw = t1 + h / 2
    xw = tw / 2

    # Bottom Flange
    A2 = b2 * t2
    y2 = t1 + h + t2 / 2
    x2 = b2 / 2

    # 3. Calculate Centroid (xc, yc)
    TotalArea = A1 + Aw + A2
    xc = (A1 * x1 + Aw * xw + A2 * x2) / TotalArea
    yc = (A1 * y1 + Aw * yw + A2 * y2) / TotalArea

    print("--- Calculated Geometric Properties ---")
    print(f"Total Area: {TotalArea:.2f} mm^2")
    print(f"Centroid Location (xc, yc): ({xc:.2f}, {yc:.2f}) mm\n")

    # 4. Calculate Moment of Inertia (Ix) about the horizontal centroidal axis
    # Using the Parallel Axis Theorem: Ix = I_local + A*d^2
    # Vertical distance from global centroid (yc) to local centroid of each part
    dy1 = yc - y1
    dyw = yc - yw
    dy2 = yc - y2

    # Moment of inertia of each part about its own centroidal axis
    Ix1_local = (b1 * t1**3) / 12
    Ixw_local = (tw * h**3) / 12
    Ix2_local = (b2 * t2**3) / 12

    # Total moment of inertia using parallel axis theorem
    Ix = (Ix1_local + A1 * dy1**2) + (Ixw_local + Aw * dyw**2) + (Ix2_local + A2 * dy2**2)
    print(f"Moment of Inertia (Ix): {Ix:.2e} mm^4\n")

    # 5. Calculate Shear Center horizontal offset from the web's centerline
    # This offset `e` makes the moment from the external shear force
    # balance the internal moment from the shear flow in the flanges.
    # The offset from the centroid is e_from_centroid.
    # e_from_centroid = (Moment from Flanges) / V
    # e_from_centroid = (1 / (2*Ix)) * (t1*b1^2*dy1^2 + t2*b2^2*dy2^2)
    
    h1_arm = abs(dy1)  # Vertical moment arm for top flange
    h2_arm = abs(dy2)  # Vertical moment arm for bottom flange

    moment_term1 = t1 * (b1**2) * (h1_arm**2)
    moment_term2 = t2 * (b2**2) * (h2_arm**2)

    e_from_centroid = (moment_term1 + moment_term2) / (2 * Ix)

    print("--- Shear Center Calculation ---")
    print("Equation for offset 'e' from centroid: e = (t1*b1^2*h1_arm^2 + t2*b2^2*h2_arm^2) / (2 * Ix)")
    print(f"t1*b1^2*h1_arm^2 = {t1:.1f} * {b1:.1f}^2 * {h1_arm:.2f}^2 = {moment_term1:.2e}")
    print(f"t2*b2^2*h2_arm^2 = {t2:.1f} * {b2:.1f}^2 * {h2_arm:.2f}^2 = {moment_term2:.2e}")
    print(f"2 * Ix = 2 * {Ix:.2e} = {2*Ix:.2e}")
    print(f"Resulting offset 'e' from centroid: {e_from_centroid:.2f} mm\n")

    # 6. Determine Final Location of Shear Center (x_sc, y_sc)
    # The vertical position of the shear center is assumed to be at the centroid's height (a common simplification).
    y_sc = yc
    # The horizontal position is offset from the web's centerline (at tw/2).
    # We calculated the offset from the centroid, so we apply it to the centroid's x-coordinate.
    x_sc = xc - e_from_centroid

    print("--- Final Locations (relative to top-left corner) ---")
    print(f"Centroid (xc, yc):          ({xc:.2f}, {yc:.2f}) mm")
    print(f"Shear Center (x_sc, y_sc):  ({x_sc:.2f}, {y_sc:.2f}) mm")
    print("\nConclusion: The shear center's x-coordinate is negative, which means it is located")
    print("to the left of the web, outside the physical cross-section, as stated in choice G.")


if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()