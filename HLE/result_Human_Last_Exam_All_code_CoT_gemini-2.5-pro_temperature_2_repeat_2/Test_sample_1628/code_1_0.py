def calculate_and_print_shear_center():
    """
    Calculates and prints the shear center location for an asymmetric channel section.

    This calculation uses a thin-walled centerline model, where:
    - h is the web height between flange centerlines.
    - b1, b2 are flange widths measured from the web centerline.
    - t_w, t_f1, t_f2 are the thicknesses.
    The shear center offset 'e' is its horizontal distance from the web's centerline.
    """
    # --- 1. Define Section Geometry (using centerline dimensions in mm) ---
    h = 150.0     # Height of web between flange centerlines
    t_w = 8.0     # Web thickness
    b1 = 50.0     # Top flange width from web centerline
    t_f1 = 10.0   # Top flange thickness
    b2 = 100.0    # Bottom flange width from web centerline
    t_f2 = 12.0   # Bottom flange thickness

    print("--- Asymmetric Channel Section Properties (Centerline Model, mm) ---")
    print(f"Web Height (h): {h}")
    print(f"Web Thickness (t_w): {t_w}")
    print(f"Top Flange Width (b1): {b1}")
    print(f"Top Flange Thickness (t_f1): {t_f1}")
    print(f"Bottom Flange Width (b2): {b2}")
    print(f"Bottom Flange Thickness (t_f2): {t_f2}")
    print("-" * 60)

    # --- 2. Calculate Vertical Centroid (y_bar) from the web's centerline ---
    A_w = h * t_w
    A_f1 = b1 * t_f1
    A_f2 = b2 * t_f2
    A_total = A_w + A_f1 + A_f2

    # Moment of area about the web's centerline (origin at web center)
    # y-coordinates are +h/2 for top flange and -h/2 for bottom flange
    moment_of_area_y = (A_f1 * (h / 2)) + (A_f2 * (-h / 2))

    y_bar = moment_of_area_y / A_total
    
    print("--- Calculation Results ---")
    print(f"Total Area (A_total): {A_total:.2f} mm^2")
    print(f"Vertical Centroid from Web Centerline (y_bar): {y_bar:.2f} mm")

    # --- 3. Calculate Moment of Inertia (I_z) about the Neutral Axis ---
    # Distances from the Neutral Axis (y_bar) to the component centroids
    d_w = 0 - y_bar  # Centroid of web is at y=0 in our chosen system
    d_f1 = (h / 2) - y_bar
    d_f2 = (-h / 2) - y_bar

    # Using Parallel Axis Theorem: I_z = I_c + A*d^2
    # For thin flanges, their own I_c (e.g., b1*t_f1**3/12) is negligible. We include it for completeness.
    I_w = (t_w * h**3 / 12) + A_w * d_w**2
    I_f1 = (b1 * t_f1**3 / 12) + A_f1 * d_f1**2
    I_f2 = (b2 * t_f2**3 / 12) + A_f2 * d_f2**2
    
    I_z = I_w + I_f1 + I_f2
    print(f"Moment of Inertia about Neutral Axis (I_z): {I_z:.2f} mm^4")

    # --- 4. Calculate Shear Center horizontal offset (e) ---
    # The formula is e = sum(Moment from flange shear forces) / (V*I_z)
    # Simplified to: e = (t_f1*h1^2*b1^2 + t_f2*h2^2*b2^2) / (2*I_z)
    # where h1 and h2 are the distances from the neutral axis to the flange centerlines.
    h1 = d_f1
    h2 = abs(d_f2)

    numerator = (t_f1 * h1**2 * b1**2) + (t_f2 * h2**2 * b2**2)
    e = numerator / (2 * I_z)

    print("-" * 60)
    print("--- Final Shear Center Location ---")
    print("The shear center is the point where shear force can be applied without causing twisting.")
    print(f"Horizontal Offset from Web Centerline (e): {e:.2f} mm")
    print(f"This location is outside the physical cross-section, on the side away from the flanges.")


# Execute the function
calculate_and_print_shear_center()