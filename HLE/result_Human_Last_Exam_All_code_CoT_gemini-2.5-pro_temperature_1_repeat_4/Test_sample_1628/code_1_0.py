import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center for an asymmetric channel section.
    The origin (0,0) of the coordinate system is set at the exterior corner
    of the bottom flange and the web's back face.
    """
    # --- Define Section Dimensions (in mm) ---
    h_web = 180.0   # Height of the web
    t_web = 10.0    # Thickness of the web
    b_f_top = 100.0  # Total width of the top flange
    t_f_top = 15.0   # Thickness of the top flange
    b_f_bottom = 60.0 # Total width of the bottom flange
    t_f_bottom = 12.0 # Thickness of the bottom flange

    print("--- Section Dimensions (mm) ---")
    print(f"Web Height (h_web): {h_web}")
    print(f"Web Thickness (t_web): {t_web}")
    print(f"Top Flange Width (b_f_top): {b_f_top}")
    print(f"Top Flange Thickness (t_f_top): {t_f_top}")
    print(f"Bottom Flange Width (b_f_bottom): {b_f_bottom}")
    print(f"Bottom Flange Thickness (t_f_bottom): {t_f_bottom}\n")

    # --- 1. Calculate Centroid (y_bar) ---
    # The section is divided into three rectangles: bottom flange, web, top flange.
    # Area and y-centroid of bottom flange
    area_bot = b_f_bottom * t_f_bottom
    y_bot = t_f_bottom / 2.0
    # Area and y-centroid of web
    area_web = h_web * t_web
    y_web = t_f_bottom + h_web / 2.0
    # Area and y-centroid of top flange
    area_top = b_f_top * t_f_top
    y_top = t_f_bottom + h_web + t_f_top / 2.0

    total_area = area_bot + area_web + area_top
    y_bar = (area_bot * y_bot + area_web * y_web + area_top * y_top) / total_area

    print("--- 1. Centroid Calculation ---")
    print(f"Centroid location from bottom edge (y_bar): {y_bar:.2f} mm\n")

    # --- 2. Calculate Moment of Inertia (I_x) about the centroidal axis ---
    # Distances from component centroids to the overall centroidal axis (y_bar)
    d_bot = y_bot - y_bar
    d_web = y_web - y_bar
    d_top = y_top - y_bar

    # Moment of inertia for each component about its own centroid
    I_c_bot = (b_f_bottom * t_f_bottom**3) / 12.0
    I_c_web = (t_web * h_web**3) / 12.0
    I_c_top = (b_f_top * t_f_top**3) / 12.0

    # Total moment of inertia using Parallel Axis Theorem: I_x = sum(I_c + A*d^2)
    I_x = (I_c_bot + area_bot * d_bot**2) + \
          (I_c_web + area_web * d_web**2) + \
          (I_c_top + area_top * d_top**2)

    print("--- 2. Moment of Inertia (I_x) Calculation ---")
    print(f"Moment of Inertia about horizontal centroidal axis (I_x): {I_x:.2f} mm^4\n")

    # --- 3. Calculate Shear Center Offset (e) ---
    # This offset 'e' is the horizontal distance from the web's centerline.
    # h1 and h2 are the distances from the neutral axis to the flange centerlines.
    h1 = abs(d_top)  # Distance to top flange centerline
    h2 = abs(d_bot)  # Distance to bottom flange centerline

    # Formula for shear center offset from web centerline:
    # e = (1 / (2*I_x)) * (t_f_top*b_f_top^2*h1^2 + t_f_bottom*b_f_bottom^2*h2^2)
    # Note: b_f here refers to the flange width used for shear flow integration, which is
    #       approximated by the total flange width for thin-walled sections.

    numerator_val = (t_f_top * b_f_top**2 * h1**2) + (t_f_bottom * b_f_bottom**2 * h2**2)
    e = numerator_val / (2 * I_x)

    print("--- 3. Shear Center Offset (e) Calculation ---")
    print("The formula for the horizontal offset 'e' from the web centerline is:")
    print("e = (t_f_top*b_f_top^2*h1^2 + t_f_bottom*b_f_bottom^2*h2^2) / (2 * I_x)\n")
    print("Plugging in the numbers:")
    print(f"e = ({t_f_top:.1f}*{b_f_top:.1f}^2*{h1:.2f}^2 + {t_f_bottom:.1f}*{b_f_bottom:.1f}^2*{h2:.2f}^2) / (2 * {I_x:.2f})")
    print(f"e = ({t_f_top * b_f_top**2 * h1**2:.2f} + {t_f_bottom * b_f_bottom**2 * h2**2:.2f}) / {2 * I_x:.2f}")
    print(f"e = {numerator_val:.2f} / {2 * I_x:.2f}")
    print(f"e = {e:.2f} mm\n")

    # --- 4. Final Shear Center Location ---
    # Vertical coordinate is the same as the centroid's y-coordinate.
    y_sc = y_bar
    # Horizontal coordinate is offset from the web centerline.
    # Web centerline x-coordinate = t_web / 2.
    # The offset 'e' is traditionally measured away from the flanges.
    x_sc = (t_web / 2.0) - e

    print("--- 4. Final Shear Center Location ---")
    print(f"The shear center is located at coordinates (x_sc, y_sc).")
    print("The origin (0,0) is at the back of the web, on the bottom edge.")
    print(f"Shear Center x-coordinate (x_sc): {x_sc:.2f} mm")
    print(f"Shear Center y-coordinate (y_sc): {y_sc:.2f} mm")
    print("\nSince the section's material is in the range x=[0, 100], a negative x_sc means")
    print("the shear center is located OUTSIDE the physical cross-section.")

calculate_shear_center_asymmetric_channel()