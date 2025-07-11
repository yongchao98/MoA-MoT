import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates and prints the shear center location for a thin-walled,
    asymmetric channel section with constant thickness.
    """
    # 1. Define Geometry (example dimensions in mm)
    h_web = 200.0  # Height of the web
    b1 = 100.0     # Width of the top flange
    b2 = 50.0      # Width of the bottom flange
    t = 10.0       # Constant thickness for web and flanges

    print("--- Section Properties ---")
    print(f"Web Height (h): {h_web} mm")
    print(f"Top Flange Width (b1): {b1} mm")
    print(f"Bottom Flange Width (b2): {b2} mm")
    print(f"Thickness (t): {t} mm\n")

    # 2. Locate Centroid (y_bar)
    # Reference axis is the centerline of the bottom flange.
    # Areas of the components
    A_top_flange = b1 * t
    A_web = h_web * t
    A_bot_flange = b2 * t
    A_total = A_top_flange + A_web + A_bot_flange

    # y_bar from the centerline of the bottom flange
    y_bar = (A_top_flange * h_web + A_web * (h_web / 2.0)) / A_total

    print("--- Centroid Calculation ---")
    print(f"Vertical distance of centroid from bottom flange centerline (y_bar): {y_bar:.2f} mm\n")

    # 3. Calculate Moment of Inertia (I_x) about the centroidal axis
    # Using the Parallel Axis Theorem: I = I_c + A*d^2
    # I_c for flanges about their strong axis ((b*t^3)/12) is neglected as it's small for thin walls.
    
    # Distance from component centroid to section centroid
    d_top = h_web - y_bar
    d_web = h_web / 2.0 - y_bar
    d_bot = -y_bar
    
    # I_x for each component
    I_x_top = A_top_flange * (d_top**2)
    I_x_web = (t * h_web**3) / 12.0 + A_web * (d_web**2)
    I_x_bot = A_bot_flange * (d_bot**2)
    
    I_x_total = I_x_top + I_x_web + I_x_bot
    
    print("--- Moment of Inertia Calculation ---")
    print(f"Moment of Inertia about neutral axis (I_x): {I_x_total:,.2f} mm^4\n")
    
    # 4. Calculate Shear Center eccentricity (e) from the web centerline
    # Using an approximate formula for thin-walled channel sections.
    # 'h' in the formula is the distance between flange centerlines, which is h_web here.
    e = (t * h_web**2) / (4 * I_x_total) * (b1**2 + b2**2)
    
    print("--- Shear Center Calculation ---")
    print(f"The shear center is located at a horizontal distance 'e' from the web centerline.")
    print(f"Eccentricity (e) = (t * h^2) / (4 * I_x) * (b1^2 + b2^2)")
    print(f"e = ({t} * {h_web}^2) / (4 * {I_x_total:,.2f}) * ({b1}^2 + {b2}^2)")
    print(f"Calculated Eccentricity (e): {e:.2f} mm\n")
    
    print("--- Conclusion ---")
    print("The shear center is located outside the cross-section, on the opposite side of the web from the flanges.")
    print(f"Its location, e = {e:.2f} mm from the web, is determined by the section's dimensions.")

# Run the calculation
calculate_shear_center_asymmetric_channel()