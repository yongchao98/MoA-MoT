import math

def calculate_shear_center():
    """
    Calculates the shear center location for an asymmetric channel section.
    
    The shear center is a point where a shear force can be applied without causing torsion.
    For a channel section, it is offset horizontally from the web.
    """
    
    # --- 1. Define Section Geometry (using centerline dimensions) ---
    # All units are in millimeters (mm)
    h = 180.0   # Height of the web (between flange centerlines)
    t_w = 8.0   # Thickness of the web
    b_1 = 100.0 # Width of the top flange
    t_1 = 12.0  # Thickness of the top flange
    b_2 = 60.0  # Width of the bottom flange
    t_2 = 10.0  # Thickness of the bottom flange

    print("--- Geometric Properties (mm) ---")
    print(f"Web Height (h): {h}")
    print(f"Web Thickness (t_w): {t_w}")
    print(f"Top Flange Width (b_1): {b_1}")
    print(f"Top Flange Thickness (t_1): {t_1}")
    print(f"Bottom Flange Width (b_2): {b_2}")
    print(f"Bottom Flange Thickness (t_2): {t_2}")
    print("-" * 35)

    # --- 2. Calculate Centroid (y_c) ---
    # Origin (y=0) is at the centerline of the bottom flange.
    # Top flange is at y=h.
    A_1 = b_1 * t_1       # Area of top flange
    y_1 = h               # y-centroid of top flange
    A_w = h * t_w         # Area of web
    y_w = h / 2.0         # y-centroid of web
    A_2 = b_2 * t_2       # Area of bottom flange
    y_2 = 0.0             # y-centroid of bottom flange
    
    total_area = A_1 + A_w + A_2
    sum_ay = (A_1 * y_1) + (A_w * y_w) + (A_2 * y_2)
    y_c = sum_ay / total_area

    print("--- Centroid Calculation ---")
    print(f"Total Area (A): {total_area:.2f} mm^2")
    print(f"Sum of (Area * y_distance) from origin: {sum_ay:.2f} mm^3")
    print(f"Vertical Centroid (y_c) from bottom flange centerline: {y_c:.2f} mm")
    print("-" * 35)

    # --- 3. Calculate Moment of Inertia (I_x) ---
    # Using Parallel Axis Theorem: I = I_local + A*d^2
    # Distance from section centroid to part centroid
    d_1 = y_1 - y_c
    d_w = y_w - y_c
    d_2 = y_2 - y_c
    
    # Moment of inertia for each part about its own centroid
    I_local_1 = (b_1 * t_1**3) / 12.0
    I_local_w = (t_w * h**3) / 12.0
    I_local_2 = (b_2 * t_2**3) / 12.0
    
    # Total moment of inertia about the section's centroidal x-axis
    I_x = (I_local_1 + A_1 * d_1**2) + \
          (I_local_w + A_w * d_w**2) + \
          (I_local_2 + A_2 * d_2**2)

    print("--- Moment of Inertia Calculation (I_x) ---")
    print(f"I_x for section: {I_x:,.2f} mm^4")
    print("-" * 35)

    # --- 4. Calculate Distances to Flange Centerlines from Centroid ---
    h_c1 = abs(d_1)  # Distance from centroid to top flange centerline
    h_c2 = abs(d_2)  # Distance from centroid to bottom flange centerline

    print("--- Distances from Centroid to Flanges ---")
    print(f"Distance to top flange centerline (h_c1): {h_c1:.2f} mm")
    print(f"Distance to bottom flange centerline (h_c2): {h_c2:.2f} mm")
    print("-" * 35)

    # --- 5. Calculate Shear Center offset (e) ---
    # Formula for horizontal offset 'e' from the web centerline:
    # e = (1 / (2 * I_x)) * (t_1*b_1^2*h_c1^2 + t_2*b_2^2*h_c2^2)
    numerator = (t_1 * b_1**2 * h_c1**2) + (t_2 * b_2**2 * h_c2**2)
    e = numerator / (2 * I_x)

    print("--- Shear Center Calculation (e) ---")
    print("Formula: e = (t_1*b_1^2*h_c1^2 + t_2*b_2^2*h_c2^2) / (2 * I_x)")
    print("Substituting values:")
    print(f"e = ({t_1:.1f}*{b_1:.1f}^2*{h_c1:.2f}^2 + {t_2:.1f}*{b_2:.1f}^2*{h_c2:.2f}^2) / (2 * {I_x:,.2f})")
    print(f"e = ({t_1 * b_1**2 * h_c1**2:,.2f} + {t_2 * b_2**2 * h_c2**2:,.2f}) / ({2 * I_x:,.2f})")
    print(f"e = {numerator:,.2f} / {2 * I_x:,.2f}")
    print("-" * 35)
    print(f"Final Result:")
    print(f"The shear center is located at a horizontal distance 'e' = {e:.2f} mm from the web's centerline.")

if __name__ == '__main__':
    calculate_shear_center()