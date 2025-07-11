import math

def calculate_shear_center_asymmetric_channel(b1, t1, b2, t2, h, tw):
    """
    Calculates the shear center location for an asymmetric channel section.

    The calculation is based on thin-walled beam theory. The origin (x=0, y=0) is at
    the intersection of the web's centerline and the section's mid-height. The shear
    center's offset 'e' is calculated along the x-axis from the web centerline.

    Args:
        b1 (float): Width of the top flange.
        t1 (float): Thickness of the top flange.
        b2 (float): Width of the bottom flange.
        t2 (float): Thickness of the bottom flange.
        h (float): Height of the web (from centerline to centerline of flanges).
        tw (float): Thickness of the web.

    Returns:
        tuple: A tuple containing (y_centroid, I_x, e), where:
               y_centroid is the centroid location along the y-axis from the mid-height.
               I_x is the moment of inertia about the centroidal x-axis.
               e is the shear center offset from the web centerline.
    """
    print("--- Input Section Properties ---")
    print(f"Top Flange Width (b1): {b1}")
    print(f"Top Flange Thickness (t1): {t1}")
    print(f"Bottom Flange Width (b2): {b2}")
    print(f"Bottom Flange Thickness (t2): {t2}")
    print(f"Web Height (h): {h}")
    print(f"Web Thickness (tw): {tw}\n")

    # 1. Calculate area of each component
    A_top = b1 * t1
    A_bot = b2 * t2
    A_web = h * tw
    A_total = A_top + A_bot + A_web

    # 2. Calculate y-coordinate of the centroid (y_c) from the mid-height
    # The top flange centerline is at y = h/2, bottom flange at y = -h/2
    y_c = (A_top * (h / 2) + A_bot * (-h / 2)) / A_total

    print("--- Intermediate Calculations ---")
    print(f"Centroid location from mid-height (y_c): {y_c:.4f}")

    # 3. Calculate distances from centroidal axis to flange centerlines
    h_top_dist = (h / 2) - y_c
    h_bot_dist = (h / 2) + y_c

    # 4. Calculate Moment of Inertia (I_x) about the centroidal x-axis
    # Using thin-walled approximation (neglecting I_c for flanges about their own axes)
    # and the parallel axis theorem.
    I_x_top = A_top * h_top_dist**2
    I_x_bot = A_bot * h_bot_dist**2
    I_x_web = (tw * h**3) / 12 + A_web * y_c**2
    I_x = I_x_top + I_x_bot + I_x_web

    print(f"Moment of Inertia about centroidal x-axis (I_x): {I_x:.4f}\n")
    
    # 5. Calculate Shear Center offset 'e' from the web centerline
    # This formula is derived from balancing the moment from shear flow in the flanges
    # against the moment of the applied shear force V*e.
    # Moment = Force_in_flange * distance_from_centroid
    # Force_in_flange = (V*t*h_dist_from_centroid * b^2) / (2*I_x)
    numerator = (t1 * (b1**2) * (h_top_dist**2)) + (t2 * (b2**2) * (h_bot_dist**2))
    e = numerator / (2 * I_x)

    print("--- Final Result ---")
    print("The formula for the shear center offset (e) from the web is:")
    print("e = (t1*b1²*h_top² + t2*b2²*h_bot²) / (2*I_x)")
    print(f"e = ({t1}*{b1}²*{h_top_dist:.2f}² + {t2}*{b2}²*{h_bot_dist:.2f}²) / (2*{I_x:.2f})")
    print(f"Shear Center offset from web (e): {e:.4f}")
    print(f"\nThe shear center is located at coordinates (x, y) = ({-e:.4f}, {y_c:.4f})")
    print("relative to the web centerline and section mid-height.")
    
    return y_c, I_x, e

# --- Example Usage ---
# Define the geometry of an asymmetric channel section (units can be mm, in, etc.)
top_flange_width = 100.0
top_flange_thickness = 10.0
bottom_flange_width = 60.0
bottom_flange_thickness = 8.0
web_height = 200.0
web_thickness = 6.0

calculate_shear_center_asymmetric_channel(
    b1=top_flange_width,
    t1=top_flange_thickness,
    b2=bottom_flange_width,
    t2=bottom_flange_thickness,
    h=web_height,
    tw=web_thickness
)