import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for an asymmetric channel section.
    The shear center is the point through which shear loads must act to produce
    bending without twisting. For a channel, this point is located outside
    the section, horizontally offset from the web.
    """
    # --- Define Section Geometry (in mm) ---
    # b1: Top flange width (projection from web)
    # t1: Top flange thickness
    # b2: Bottom flange width (projection from web)
    # t2: Bottom flange thickness
    # h_total: Total height of the section
    # tw: Web thickness
    
    b1 = 100.0
    t1 = 12.0
    b2 = 80.0
    t2 = 10.0
    h_total = 200.0
    tw = 8.0

    print("--- Section Dimensions (mm) ---")
    print(f"Top Flange Width (b1): {b1}")
    print(f"Top Flange Thickness (t1): {t1}")
    print(f"Bottom Flange Width (b2): {b2}")
    print(f"Bottom Flange Thickness (t2): {t2}")
    print(f"Total Height (h_total): {h_total}")
    print(f"Web Thickness (tw): {tw}\n")

    # 1. Calculate Centroid Location (y_bar) from the bottom outer edge
    h_web = h_total - t1 - t2

    # Areas of the three rectangular parts
    A_top = b1 * t1
    A_web = h_web * tw
    A_bot = b2 * t2
    A_total = A_top + A_web + A_bot

    # Local y-centroids of each part (from bottom outer edge)
    y_top_local = h_total - t1 / 2.0
    y_web_local = t2 + h_web / 2.0
    y_bot_local = t2 / 2.0

    # Overall centroid y_bar
    y_bar = (A_top * y_top_local + A_web * y_web_local + A_bot * y_bot_local) / A_total

    print("--- Intermediate Calculations ---")
    print(f"Total Area (A_total): {A_total:.2f} mm^2")
    print(f"Centroid location (y_bar) from bottom edge: {y_bar:.2f} mm\n")

    # 2. Calculate Moment of Inertia (I_x) about the horizontal centroidal axis
    # Using the Parallel Axis Theorem: I = I_local + A*d^2
    
    # Moment of inertia for top flange
    I_x_top = (b1 * t1**3 / 12.0) + A_top * (y_top_local - y_bar)**2
    
    # Moment of inertia for web
    I_x_web = (tw * h_web**3 / 12.0) + A_web * (y_web_local - y_bar)**2
    
    # Moment of inertia for bottom flange
    I_x_bot = (b2 * t2**3 / 12.0) + A_bot * (y_bot_local - y_bar)**2

    I_x_total = I_x_top + I_x_web + I_x_bot
    
    print(f"Moment of Inertia (I_x_total): {I_x_total:.2f} mm^4\n")

    # 3. Calculate Shear Center Eccentricity (e)
    # This is the horizontal distance from the web's centerline.

    # Distances from the section's centroid to the centerline of each flange
    d1 = y_top_local - y_bar  # Distance to top flange centerline
    d2 = y_bar - y_bot_local  # Distance to bottom flange centerline

    # Formula for shear center eccentricity 'e'
    # e = (1 / (2*I_x)) * (t1*d1^2*b1^2 + t2*d2^2*b2^2)
    numerator = t1 * (d1**2) * (b1**2) + t2 * (d2**2) * (b2**2)
    denominator = 2 * I_x_total
    e = numerator / denominator

    print("--- Final Result ---")
    print(f"The shear center is located at a horizontal distance 'e' from the centerline of the web.")
    print(f"Calculated eccentricity e = ({t1:.1f} * {d1:.2f}^2 * {b1:.1f}^2 + {t2:.1f} * {d2:.2f}^2 * {b2:.1f}^2) / (2 * {I_x_total:.2f})")
    print(f"e = {e:.2f} mm")

if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()