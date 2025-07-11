import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for an asymmetric channel section.
    The shear center is the point where a vertical shear force can be applied
    without causing torsion. For a channel, this location 'e' is measured
    horizontally from the centerline of the web.
    """
    # --- Define Section Dimensions (using a sample asymmetric channel) ---
    # All units are in millimeters (mm)
    H_overall = 200.0  # Overall height
    B1_overall = 80.0   # Overall width of the top flange
    B2_overall = 120.0  # Overall width of the bottom flange
    t_w = 8.0         # Thickness of the web
    t_f1 = 12.0       # Thickness of the top flange
    t_f2 = 12.0       # Thickness of the bottom flange

    print("--- Input Section Dimensions (mm) ---")
    print(f"Overall Height (H): {H_overall}")
    print(f"Top Flange Width (B1): {B1_overall}")
    print(f"Bottom Flange Width (B2): {B2_overall}")
    print(f"Web Thickness (tw): {t_w}")
    print(f"Top Flange Thickness (tf1): {t_f1}")
    print(f"Bottom Flange Thickness (tf2): {t_f2}")
    print("-" * 35)

    # --- Step 1: Calculate Component Properties ---
    # We will use an origin at the bottom-left outer corner.
    # Component 1: Top Flange
    A1 = B1_overall * t_f1
    y1 = H_overall - t_f1 / 2.0

    # Component 2: Web
    h_web = H_overall - t_f1 - t_f2
    A2 = h_web * t_w
    y2 = h_web / 2.0 + t_f2

    # Component 3: Bottom Flange
    A3 = B2_overall * t_f2
    y3 = t_f2 / 2.0

    # --- Step 2: Locate the Centroid (y-coordinate) ---
    A_total = A1 + A2 + A3
    y_c = (A1 * y1 + A2 * y2 + A3 * y3) / A_total

    print("--- Centroid Calculation ---")
    print(f"Total Area (A_total): {A_total:.2f} mm^2")
    print(f"Centroid location from bottom (y_c): {y_c:.2f} mm")
    print("-" * 35)

    # --- Step 3: Calculate Moment of Inertia (I_x) about the centroidal axis ---
    # Using the Parallel Axis Theorem: I = I_local + A*d^2
    # For Top Flange
    I_x1_local = (B1_overall * t_f1**3) / 12.0
    d1 = y1 - y_c
    I_x1 = I_x1_local + A1 * d1**2

    # For Web
    I_x2_local = (t_w * h_web**3) / 12.0
    d2 = y2 - y_c
    I_x2 = I_x2_local + A2 * d2**2

    # For Bottom Flange
    I_x3_local = (B2_overall * t_f2**3) / 12.0
    d3 = y3 - y_c
    I_x3 = I_x3_local + A3 * d3**2

    I_x_total = I_x1 + I_x2 + I_x3

    print("--- Moment of Inertia Calculation ---")
    print(f"Moment of Inertia about centroidal x-axis (Ix): {I_x_total:.2f} mm^4")
    print("-" * 35)

    # --- Step 4: Calculate Shear Center Offset (e) from the web centerline ---
    # Flange projections from web centerline
    b1_proj = B1_overall - t_w / 2.0
    b2_proj = B2_overall - t_w / 2.0

    # Distances from centroid to flange centerlines
    y_top = d1
    y_bot = d3 # Note: y_bot is negative since the centroid is above it.

    # Numerator of the formula for 'e'
    numerator = (t_f1 * b1_proj**2 * y_top**2) + (t_f2 * b2_proj**2 * y_bot**2)
    
    # The formula requires dividing by (2 * I_x)
    e = numerator / (2 * I_x_total)

    print("--- Shear Center Calculation ---")
    print(f"Formula: e = (t_f1*b1_proj^2*y_top^2 + t_f2*b2_proj^2*y_bot^2) / (2 * Ix)")
    print(f"Calculation: e = ({t_f1}*{b1_proj:.2f}^2*{y_top:.2f}^2 + {t_f2}*{b2_proj:.2f}^2*{y_bot:.2f}^2) / (2 * {I_x_total:.2f})")
    print(f"Result: The shear center is located at a distance 'e' = {e:.2f} mm")
    print("This distance is measured horizontally from the centerline of the web, away from the flanges.")
    print("-" * 35)


if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()