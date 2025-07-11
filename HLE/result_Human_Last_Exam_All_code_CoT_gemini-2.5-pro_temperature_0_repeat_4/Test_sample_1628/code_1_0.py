import math

def calculate_asymmetric_channel_shear_center():
    """
    Calculates and explains the shear center location for an asymmetric channel section.

    The shear center is the point through which a shear force can be applied
    without causing torsion. For a channel section subjected to a vertical
    shear force, the primary interest is the horizontal offset 'e' of the
    shear center from the web.

    This script calculates:
    1. The vertical position of the centroid (Neutral Axis).
    2. The moment of inertia (I_z) about the horizontal centroidal axis.
    3. The horizontal distance 'e' of the shear center from the web's centerline.
    """

    # --- Define Section Dimensions (example values) ---
    # All units are in millimeters (mm)
    b1 = 60.0   # Top flange width
    t1 = 10.0   # Top flange thickness
    b2 = 80.0   # Bottom flange width
    t2 = 12.0   # Bottom flange thickness
    h = 150.0   # Web height (clear height between flanges)
    tw = 8.0    # Web thickness

    print("--- Section Dimensions ---")
    print(f"Top Flange Width (b1): {b1} mm")
    print(f"Top Flange Thickness (t1): {t1} mm")
    print(f"Bottom Flange Width (b2): {b2} mm")
    print(f"Bottom Flange Thickness (t2): {t2} mm")
    print(f"Web Height (h): {h} mm")
    print(f"Web Thickness (tw): {tw} mm")
    print("-" * 35 + "\n")

    # --- 1. Calculate Centroid (Neutral Axis location) ---
    # Origin (y=0) is set at the bottom face of the bottom flange.
    
    # y-coordinate of the centroid of each component
    y_comp_top_flange = t2 + h + t1 / 2.0
    y_comp_web = t2 + h / 2.0
    y_comp_bottom_flange = t2 / 2.0
    
    # Area of each component
    A_top_flange = b1 * t1
    A_web = h * tw
    A_bottom_flange = b2 * t2
    A_total = A_top_flange + A_web + A_bottom_flange
    
    # Vertical centroid (y_c) from the bottom face
    y_c = ((A_top_flange * y_comp_top_flange) + (A_web * y_comp_web) + (A_bottom_flange * y_comp_bottom_flange)) / A_total

    print("--- Centroid Calculation ---")
    print(f"Vertical centroid (y_c) from bottom face: {y_c:.2f} mm")
    print("-" * 35 + "\n")

    # --- 2. Calculate Moment of Inertia (I_z) ---
    # Using the Parallel Axis Theorem: I_z = I_zc + A*d^2
    
    # Distances from component centroids to the section's neutral axis
    d_top_flange = y_comp_top_flange - y_c
    d_web = y_comp_web - y_c
    d_bottom_flange = y_comp_bottom_flange - y_c

    # I_z for each component about the section's neutral axis
    Iz_top_flange = (b1 * t1**3 / 12.0) + A_top_flange * d_top_flange**2
    Iz_web = (tw * h**3 / 12.0) + A_web * d_web**2
    Iz_bottom_flange = (b2 * t2**3 / 12.0) + A_bottom_flange * d_bottom_flange**2

    # Total Moment of Inertia
    I_z = Iz_top_flange + Iz_web + Iz_bottom_flange

    print("--- Moment of Inertia Calculation ---")
    print(f"Moment of Inertia (I_z): {I_z:,.2f} mm^4")
    print("-" * 35 + "\n")

    # --- 3. Calculate Shear Center horizontal offset 'e' ---
    # The formula is derived from balancing the external moment (V*e) with the
    # internal moment from the shear flow in the flanges.
    # e = (1 / (2*I_z)) * [ t1 * b1^2 * d_top^2 + t2 * b2^2 * d_bottom^2 ]
    # where d_top and d_bottom are the absolute distances from the neutral axis
    # to the respective flange centerlines.
    
    e = (1 / (2 * I_z)) * (t1 * b1**2 * d_top_flange**2 + t2 * b2**2 * abs(d_bottom_flange)**2)

    print("--- Shear Center Calculation ---")
    print(f"Distance from top flange centerline to NA: {d_top_flange:.2f} mm")
    print(f"Distance from bottom flange centerline to NA: {abs(d_bottom_flange):.2f} mm")
    print(f"Horizontal offset 'e' of Shear Center from web centerline: {e:.2f} mm")
    print("\n" + "="*35)
    print("             CONCLUSION")
    print("="*35)
    print("The shear center is located OUTSIDE the cross-section.")
    print(f"It is offset horizontally from the web by a calculated distance e = {e:.2f} mm.")
    print("This location ensures that a vertical shear force applied through it will not cause twisting.")
    print("This confirms the location is determined by the channel's specific dimensions.")
    print("="*35)

# Execute the function
calculate_asymmetric_channel_shear_center()