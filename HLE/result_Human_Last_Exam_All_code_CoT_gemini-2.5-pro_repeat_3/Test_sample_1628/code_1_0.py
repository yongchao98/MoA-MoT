import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for an asymmetric channel section.
    
    The shear center is the point where a shear force can be applied without
    causing torsion. For a channel, it is located at a horizontal distance 'e'
    from the web.
    """
    
    # --- 1. Define Geometry (example dimensions in inches) ---
    H = 10.0   # Total height of the section
    b1 = 3.0   # Top flange width (from web face)
    b2 = 4.0   # Bottom flange width (from web face)
    tw = 0.5   # Web thickness
    t1 = 0.5   # Top flange thickness
    t2 = 0.6   # Bottom flange thickness

    print("1. Given Section Dimensions (in):")
    print(f"   - Total Height (H): {H}")
    print(f"   - Top Flange Width (b1): {b1}")
    print(f"   - Bottom Flange Width (b2): {b2}")
    print(f"   - Web Thickness (tw): {tw}")
    print(f"   - Top Flange Thickness (t1): {t1}")
    print(f"   - Bottom Flange Thickness (t2): {t2}\n")

    # Define the three rectangular components
    # Origin (0,0) is at the bottom-left exterior corner
    # Part 1: Top Flange
    A1 = (b1 + tw) * t1 
    y1 = H - t1 / 2.0
    # Part 2: Web
    A_web = (H - t1 - t2) * tw
    y_web = (H - t1 - t2) / 2.0 + t2
    # Part 3: Bottom Flange
    A2 = (b2 + tw) * t2
    y2 = t2 / 2.0
    
    # Total Area
    A_total = A1 + A_web + A2

    # --- 2. Locate the Centroid (Cy) ---
    # Calculate the vertical position of the centroid from the bottom edge
    Cy = (A1 * y1 + A_web * y_web + A2 * y2) / A_total

    print("2. Centroid Calculation:")
    print(f"   - Total Area (A_total): {A_total:.4f} in^2")
    print(f"   - Vertical Centroid (Cy) from bottom: {Cy:.4f} in\n")

    # --- 3. Calculate Moment of Inertia (Ix) about the centroidal axis ---
    # Use Parallel Axis Theorem: I = I_c + A*d^2
    # For a rectangle, I_c about its own centroid is (b*h^3)/12
    
    # Top Flange
    Ic1 = ((b1 + tw) * t1**3) / 12.0
    d1_centroid = y1 - Cy
    Ix1 = Ic1 + A1 * d1_centroid**2

    # Web
    Ic_web = (tw * (H - t1 - t2)**3) / 12.0
    d_web_centroid = y_web - Cy
    Ix_web = Ic_web + A_web * d_web_centroid**2

    # Bottom Flange
    Ic2 = ((b2 + tw) * t2**3) / 12.0
    d2_centroid = y2 - Cy
    Ix2 = Ic2 + A2 * d2_centroid**2

    I_x = Ix1 + Ix_web + Ix2

    print("3. Moment of Inertia Calculation:")
    print(f"   - Moment of Inertia (I_x): {I_x:.4f} in^4\n")

    # --- 4. Calculate Shear Center Offset (e) ---
    # The horizontal offset 'e' from the web's centerline.
    # We use a simplified model considering flanges as rectangles extending from web centerline.
    # Widths are measured from web centerline for this formula.
    flange_b1 = b1 + tw / 2.0
    flange_b2 = b2 + tw / 2.0
    # d1 and d2 are distances from the section's centroid to the flange's centerline
    d1 = H - t1 / 2.0 - Cy  # Distance for top flange
    d2 = Cy - t2 / 2.0      # Distance for bottom flange
    
    # Formula for 'e' based on moment balance of shear flow
    numerator = (t1 * flange_b1**2 * d1**2) + (t2 * flange_b2**2 * d2**2)
    e = numerator / (2 * I_x)

    print("4. Shear Center Calculation:")
    print(f"   - Distance from centroid to top flange centerline (d1): {d1:.4f} in")
    print(f"   - Distance from centroid to bottom flange centerline (d2): {d2:.4f} in")
    print(f"   - Horizontal offset 'e' from web centerline: {e:.4f} in\n")
    
    print("Conclusion:")
    print("The shear center is located outside the section,")
    print(f"at a horizontal distance e = {e:.4f} inches from the web's centerline,")
    print(f"and at a vertical height of Cy = {Cy:.4f} inches from the bottom.")

if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()