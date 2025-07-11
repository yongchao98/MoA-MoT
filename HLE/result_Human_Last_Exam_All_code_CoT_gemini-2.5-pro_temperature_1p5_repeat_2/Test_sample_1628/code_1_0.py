import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for an asymmetric channel section.

    The shear center is the point through which shear loads must act to produce
    bending without twisting. For a channel section, it is located outside the
    section at an eccentricity 'e' from the web.

    This script calculates 'e' based on the section's geometric properties.
    """
    # --- 1. Define Section Properties (using centerline dimensions) ---
    # All units are in millimeters (mm)
    b1 = 60.0   # Width of the top flange
    t1 = 8.0    # Thickness of the top flange
    b2 = 100.0  # Width of the bottom flange
    t2 = 10.0   # Thickness of the bottom flange
    h = 200.0   # Height of the web (between flange centerlines)
    tw = 6.0    # Thickness of the web

    print("--- Section Properties (mm) ---")
    print(f"Top Flange Width (b1): {b1}, Thickness (t1): {t1}")
    print(f"Bottom Flange Width (b2): {b2}, Thickness (t2): {t2}")
    print(f"Web Height (h): {h}, Thickness (tw): {tw}\n")

    # --- 2. Calculate Centroid Location (y_bar) ---
    # Measured from the centerline of the bottom flange
    A1 = b1 * t1  # Area of top flange
    A2 = b2 * t2  # Area of bottom flange
    Aw = h * tw   # Area of web

    # Sum of (Area * y_distance_from_ref) / Sum of Areas
    y_bar = (A1 * h + Aw * (h / 2)) / (A1 + A2 + Aw)

    print("--- Intermediate Calculations ---")
    print(f"Centroid location (y_bar) from bottom flange centerline: {y_bar:.2f} mm\n")

    # --- 3. Calculate Moment of Inertia (I_x) about the horizontal centroidal axis ---
    # Using the Parallel Axis Theorem: I = I_local + A*d^2
    # For thin rectangles, I_local (e.g., b*t^3/12) is small and often ignored, but we will include it for accuracy.

    # Distances of component centroids from the overall section centroid
    d1 = h - y_bar  # Distance from centroid to top flange centerline
    d2 = -y_bar     # Distance from centroid to bottom flange centerline
    dw = h / 2 - y_bar # Distance from centroid to web centerline

    I_x1 = (b1 * t1**3 / 12) + A1 * d1**2  # I_x for top flange
    I_x2 = (b2 * t2**3 / 12) + A2 * d2**2  # I_x for bottom flange
    I_xw = (tw * h**3 / 12) + Aw * dw**2  # I_x for web

    I_x = I_x1 + I_x2 + I_xw

    print(f"Moment of Inertia about horizontal centroidal axis (I_x): {I_x:,.2f} mm^4\n")

    # --- 4. Calculate Shear Center Eccentricity (e) ---
    # e is the distance from the web's centerline
    # Formula derived from balancing external and internal moments
    
    h1 = h - y_bar # Vertical distance from centroid to top flange centerline
    h2 = y_bar     # Vertical distance from centroid to bottom flange centerline

    # Moment contribution from top flange shear flow
    moment_top = t1 * (h1**2) * (b1**2)

    # Moment contribution from bottom flange shear flow
    moment_bottom = t2 * (h2**2) * (b2**2)
    
    e = (1 / (2 * I_x)) * (moment_top + moment_bottom)

    print("--- Final Result ---")
    print(f"The shear center is located at a distance 'e' from the centerline of the web.")
    print(f"Equation: e = (1 / (2 * I_x)) * (t1 * h1^2 * b1^2 + t2 * h2^2 * b2^2)")
    print(f"e = (1 / (2 * {I_x:,.2f})) * ({t1:.1f} * {h1:.2f}^2 * {b1:.1f}^2 + {t2:.1f} * {h2:.2f}^2 * {b2:.1f}^2)")
    print(f"Calculated Shear Center Eccentricity (e): {e:.2f} mm")
    print("\nThis location is outside the physical cross-section, offset from the web.")

calculate_shear_center_asymmetric_channel()