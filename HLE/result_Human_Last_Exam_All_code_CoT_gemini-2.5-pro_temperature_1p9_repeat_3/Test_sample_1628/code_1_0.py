import math

def calculate_asymmetric_channel_shear_center():
    """
    Calculates the shear center location for an asymmetric (unequal flange)
    channel section.

    The shear center is the point through which a shear force can be applied
    without causing torsion (twisting). For a channel section, it is located
    horizontally outside the section, offset from the web.

    This calculation assumes a thin-walled section.
    """
    # --- Define Section Dimensions (User can modify these values) ---
    # Note: Flange widths are defined as the length projecting from the web face.
    b1 = 60  # Top flange projecting width (mm)
    t1 = 10  # Top flange thickness (mm)
    b2 = 40  # Bottom flange projecting width (mm)
    t2 = 8   # Bottom flange thickness (mm)
    h_web = 150 # Height of the web (mm)
    tw = 6      # Web thickness (mm)

    print("--- Input Section Dimensions ---")
    print(f"Top flange width (b1): {b1} mm")
    print(f"Top flange thickness (t1): {t1} mm")
    print(f"Bottom flange width (b2): {b2} mm")
    print(f"Bottom flange thickness (t2): {t2} mm")
    print(f"Web height (h_web): {h_web} mm")
    print(f"Web thickness (tw): {tw} mm\n")

    # --- Step 1: Define component properties relative to a common origin ---
    # Origin (y=0) is set at the centerline of the bottom flange.
    
    # Top Flange
    A1 = b1 * t1
    y1 = h_web + t1 / 2  # Centroid of top flange

    # Web
    A_web = h_web * tw
    y_web = h_web / 2   # Centroid of web (from bottom flange centerline)

    # Bottom Flange
    A2 = b2 * t2
    y2 = 0              # Centroid of bottom flange is at the origin

    # --- Step 2: Calculate the vertical centroid (Neutral Axis) of the entire section ---
    total_area = A1 + A_web + A2
    y_bar = (A1 * y1 + A_web * y_web + A2 * y2) / total_area

    print("--- Intermediate Calculations ---")
    print(f"Vertical centroid location (y_bar) from bottom flange centerline: {y_bar:.2f} mm")

    # --- Step 3: Calculate Moment of Inertia (I_x) about the centroidal axis ---
    # Using the Parallel Axis Theorem: I_x = Σ(I_c + A*d^2)
    
    # Top Flange Contribution
    Ic1 = (b1 * t1**3) / 12
    d1 = y1 - y_bar
    I_x1 = Ic1 + A1 * d1**2

    # Web Contribution
    Ic_web = (tw * h_web**3) / 12
    d_web = y_web - y_bar
    I_x_web = Ic_web + A_web * d_web**2
    
    # Bottom Flange Contribution
    Ic2 = (b2 * t2**3) / 12
    d2 = y2 - y_bar
    I_x2 = Ic2 + A2 * d2**2

    I_x = I_x1 + I_x_web + I_x2
    print(f"Moment of Inertia about neutral x-axis (Ix): {I_x:.2e} mm^4")

    # --- Step 4 & 5: Calculate shear center offset 'e' from the web centerline ---
    # This is derived by balancing the external moment (V*e) with the internal
    # moment from the shear flow in the flanges.
    # Formula derived from first principles: e = (1/(2*Ix)) * (A1*b1*d1^2 + A2*b2*d2^2)
    # A1 = Area of top flange, b1 = width of top flange, d1 = dist of top flange centroid from NA
    # A2 = Area of bottom flange, b2 = width of bottom flange, d2 = dist of bottom flange centroid from NA
    
    # Note: We use the squared distances (d1^2, d2^2) as they appear in the derivation
    # for the moment created by the shear forces in the flanges.
    
    numerator = (A1 * b1 * d1**2) + (A2 * b2 * d2**2)
    e = numerator / (2 * I_x)

    print("\n--- Final Result ---")
    print("The shear center is a property of the cross-section.")
    print("Its location is calculated to be offset from the web centerline.")
    print(f"Equation for the offset 'e': e = (A_top*b_top*d_top^2 + A_bot*b_bot*d_bot^2) / (2 * Ix)")
    print(f"e = ({A1:.0f}*{b1:.0f}*{d1:.2f}^2 + {A2:.0f}*{b2:.0f}*{d2:.2f}^2) / (2 * {I_x:.2e})")
    print(f"Calculated shear center horizontal offset (e) from the web centerline: {e:.2f} mm")


if __name__ == '__main__':
    calculate_asymmetric_channel_shear_center()