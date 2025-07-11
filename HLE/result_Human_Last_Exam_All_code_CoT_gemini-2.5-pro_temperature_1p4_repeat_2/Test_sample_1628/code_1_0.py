import math

def calculate_shear_center():
    """
    Calculates the shear center location for an asymmetric channel section.
    The shear center is a property of the cross-section where a shear
    force can be applied without causing torsion. For a channel section,
    it is located outside the section, offset horizontally from the web.
    """

    # --- 1. Define Section Dimensions (in mm) ---
    # Top flange width
    b1 = 100.0
    # Bottom flange width
    b2 = 50.0
    # Web height (clear distance between flanges)
    h_w = 200.0
    # Flange thickness (assumed constant for both flanges)
    t_f = 12.0
    # Web thickness
    t_w = 8.0

    print("--- Asymmetric Channel Section Properties ---")
    print(f"Top Flange Width (b1): {b1:.2f} mm")
    print(f"Bottom Flange Width (b2): {b2:.2f} mm")
    print(f"Web Height (h_w): {h_w:.2f} mm")
    print(f"Flange Thickness (t_f): {t_f:.2f} mm")
    print(f"Web Thickness (t_w): {t_w:.2f} mm\n")

    # --- 2. Calculate Centroid (y_c) from bottom edge ---
    # Define the three parts: top flange, web, bottom flange
    parts = [
        {'area': b1 * t_f, 'y_bar': t_f + h_w + t_f / 2.0},  # Top flange
        {'area': h_w * t_w, 'y_bar': t_f + h_w / 2.0},      # Web
        {'area': b2 * t_f, 'y_bar': t_f / 2.0}               # Bottom flange
    ]

    total_area = sum(p['area'] for p in parts)
    sum_ay = sum(p['area'] * p['y_bar'] for p in parts)
    y_c = sum_ay / total_area

    print("--- Intermediate Calculations ---")
    print(f"Total Area (A): {total_area:.2f} mm^2")
    print(f"Centroid Location from bottom (y_c): {y_c:.2f} mm\n")

    # --- 3. Calculate Moment of Inertia (I_x) about Centroidal Axis ---
    # Using the Parallel Axis Theorem: I_x = sum(I_c + A*d^2)
    # I_c for a rectangle = (b*h^3)/12
    I_x = 0
    # Top Flange
    I_c_tf = (b1 * t_f**3) / 12.0
    d_tf = parts[0]['y_bar'] - y_c
    I_x += I_c_tf + parts[0]['area'] * d_tf**2
    # Web
    I_c_w = (t_w * h_w**3) / 12.0
    d_w = parts[1]['y_bar'] - y_c
    I_x += I_c_w + parts[1]['area'] * d_w**2
    # Bottom Flange
    I_c_bf = (b2 * t_f**3) / 12.0
    d_bf = parts[2]['y_bar'] - y_c
    I_x += I_c_bf + parts[2]['area'] * d_bf**2
    
    print(f"Moment of Inertia about x-axis (I_x): {I_x:.2f} mm^4")

    # --- 4. Calculate Distances from Centroid to Flange Centerlines ---
    y1_c = abs(d_tf)  # Distance from NA to top flange centerline
    y2_c = abs(d_bf)  # Distance from NA to bottom flange centerline
    print(f"Distance from NA to top flange centerline (y1_c): {y1_c:.2f} mm")
    print(f"Distance from NA to bottom flange centerline (y2_c): {y2_c:.2f} mm\n")
    
    # --- 5. Calculate Shear Center Eccentricity (e) ---
    # Formula for horizontal offset from the web's centerline
    # e = (t_f / (2 * I_x)) * (y1_c^2 * b1^2 + y2_c^2 * b2^2)
    e = (t_f / (2 * I_x)) * (y1_c**2 * b1**2 + y2_c**2 * b2**2)
    
    # --- 6. Print the Final Equation and Result ---
    print("--- Final Calculation for Shear Center Eccentricity (e) ---")
    print("The formula for the shear center offset 'e' from the web is:")
    print("e = (t_f / (2 * I_x)) * (y1_c^2 * b1^2 + y2_c^2 * b2^2)\n")

    print("Plugging in the calculated values:")
    final_eq_str = (
        f"e = ({t_f:.2f} / (2 * {I_x:.2f})) * "
        f"({y1_c:.2f}^2 * {b1:.2f}^2 + {y2_c:.2f}^2 * {b2:.2f}^2)"
    )
    print(final_eq_str)

    print(f"\nResult:")
    print(f"The shear center is located at a distance 'e' = {e:.2f} mm from the web.")

if __name__ == '__main__':
    calculate_shear_center()