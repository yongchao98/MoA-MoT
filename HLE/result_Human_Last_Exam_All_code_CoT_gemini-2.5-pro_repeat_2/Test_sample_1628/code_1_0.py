import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center for a thin-walled asymmetric channel section.
    The shear center is the point through which shear loads must act to produce
    bending without twisting. For a channel section, it is located at a distance 'e'
    from the web.
    """
    # 1. Define Geometry (using example dimensions in mm)
    # You can change these values to match your specific section
    h_web = 173.0  # Clear height of the web
    tw = 10.0     # Thickness of the web
    b1 = 100.0    # Width of the top flange
    t1 = 15.0     # Thickness of the top flange
    b2 = 50.0     # Width of the bottom flange
    t2 = 12.0     # Thickness of the bottom flange

    print("--- Input Section Dimensions (mm) ---")
    print(f"Web Height (h_web): {h_web}")
    print(f"Web Thickness (tw): {tw}")
    print(f"Top Flange Width (b1): {b1}")
    print(f"Top Flange Thickness (t1): {t1}")
    print(f"Bottom Flange Width (b2): {b2}")
    print(f"Bottom Flange Thickness (t2): {t2}\n")

    # 2. Locate the Centroid (y_bar)
    # We will set the origin (y=0) at the centerline of the web's bottom edge.
    # The x=0 axis is the centerline of the web.
    
    # Properties of each component part (top flange, web, bottom flange)
    A1 = b1 * t1  # Area of top flange
    y1_c = h_web + t1 / 2.0  # y-centroid of top flange

    A2 = h_web * tw  # Area of web
    y2_c = h_web / 2.0  # y-centroid of web

    A3 = b2 * t2  # Area of bottom flange
    y3_c = -t2 / 2.0  # y-centroid of bottom flange

    A_total = A1 + A2 + A3
    
    # Overall centroid y_bar relative to the web's bottom edge
    y_bar = (A1 * y1_c + A2 * y2_c + A3 * y3_c) / A_total

    print("--- Calculated Section Properties ---")
    print(f"Total Area (A): {A_total:.2f} mm^2")
    print(f"Centroid location (y_bar) from web bottom edge: {y_bar:.2f} mm\n")

    # 3. Calculate Moment of Inertia (I_x) about the centroidal axis
    # Using the Parallel Axis Theorem: I = I_c + A * d^2
    
    # For top flange
    I_c1 = (b1 * t1**3) / 12.0
    d1 = y1_c - y_bar
    I_x1 = I_c1 + A1 * d1**2

    # For web
    I_c2 = (tw * h_web**3) / 12.0
    d2 = y2_c - y_bar
    I_x2 = I_c2 + A2 * d2**2

    # For bottom flange
    I_c3 = (b2 * t2**3) / 12.0
    d3 = y3_c - y_bar
    I_x3 = I_c3 + A3 * d3**2

    I_x = I_x1 + I_x2 + I_x3
    
    print(f"Moment of Inertia about horizontal centroidal axis (Ix): {I_x:.2f} mm^4\n")

    # 4. Calculate the Shear Center offset (e) from the web centerline
    # The formula is derived from balancing the internal torque from flange shear
    # flows with the external torque V*e.
    # e = (h_flanges / (4 * I_x)) * [t1*b1^2*(y1_c - y_bar) + t2*b2^2*(y_bar - y3_c)]
    # where h_flanges is the distance between the centerlines of the flanges.
    
    h_flanges = y1_c - y3_c
    
    # Note: (y1_c - y_bar) is d1, and (y_bar - y3_c) is -d3.
    term_top_flange = t1 * b1**2 * d1
    term_bottom_flange = t2 * b2**2 * (-d3)
    
    e = (h_flanges / (4 * I_x)) * (term_top_flange + term_bottom_flange)

    print("--- Final Shear Center Location ---")
    print(f"The shear center is offset horizontally from the web's centerline.")
    print(f"Equation for offset 'e': e = (h_flanges / (4 * Ix)) * [t1*b1^2*(y1_c - y_bar) + t2*b2^2*(y_bar - y3_c)]")
    print(f"e = ({h_flanges:.2f} / (4 * {I_x:.2f})) * [{t1}*{b1}^2*({d1:.2f}) + {t2}*{b2}^2*({-d3:.2f})]")
    print(f"Offset distance (e): {e:.2f} mm")
    print(f"The shear center coordinates (e, y_bar_from_origin) are ({e:.2f}, {y_bar:.2f})")
    print("(Coordinates are relative to the web's centerline and its bottom edge)")


if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()