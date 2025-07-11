import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates and explains the shear center location for an asymmetric channel section.
    
    The shear center is the point through which an applied shear force results in bending
    without twisting. For an asymmetric channel, this point lies outside the section,
    on the side of the web, at a calculated distance 'e' from the web.
    
    This script demonstrates the calculation for a sample section.
    """
    
    # --- 1. Define Section Dimensions (in mm) ---
    # b1: Width of the top flange
    # b2: Width of the bottom flange
    # h: Height of the web (clear distance between flanges)
    # tf: Thickness of the flanges
    # tw: Thickness of the web
    b1 = 100.0
    b2 = 60.0
    h = 200.0
    tf = 10.0
    tw = 8.0

    print("--- Asymmetric Channel Section Properties ---")
    print(f"Top Flange Width (b1):    {b1} mm")
    print(f"Bottom Flange Width (b2): {b2} mm")
    print(f"Web Height (h):           {h} mm")
    print(f"Flange Thickness (tf):    {tf} mm")
    print(f"Web Thickness (tw):       {tw} mm\n")

    # --- 2. Calculate Centroid (y_c from bottom edge of bottom flange) ---
    # Area of each component
    area_flange1 = b1 * tf
    area_flange2 = b2 * tf
    area_web = h * tw
    total_area = area_flange1 + area_flange2 + area_web
    
    # Moment of area about the bottom edge
    moment_area_flange1 = area_flange1 * (tf / 2 + h + tf)
    moment_area_flange2 = area_flange2 * (tf / 2)
    moment_area_web = area_web * (h / 2 + tf)
    
    y_c = (moment_area_flange1 + moment_area_flange2 + moment_area_web) / total_area
    
    # --- 3. Calculate Moment of Inertia (I_x) about the Centroidal Axis ---
    # Using the Parallel Axis Theorem: I = I_local + A*d^2
    # Distance 'd' is the distance from the component's centroid to the section's centroid
    
    # I_x for top flange
    d1 = (tf / 2 + h + tf) - y_c
    I_x1 = (b1 * tf**3 / 12) + (area_flange1 * d1**2)

    # I_x for bottom flange
    d2 = y_c - (tf / 2)
    I_x2 = (b2 * tf**3 / 12) + (area_flange2 * d2**2)
    
    # I_x for web
    d_web = (h / 2 + tf) - y_c
    I_x_web = (tw * h**3 / 12) + (area_web * d_web**2)

    I_x_total = I_x1 + I_x2 + I_x_web

    print("--- Calculated Geometric Properties ---")
    print(f"Centroid location (y_c from bottom): {y_c:.2f} mm")
    print(f"Moment of Inertia (I_x):             {I_x_total:.2f} mm^4\n")

    # --- 4. Calculate Shear Center offset 'e' from the web centerline ---
    # This formula is derived from balancing the external shear moment (V*e)
    # with the internal moment from the shear flow in the flanges.
    # We use a thin-walled approximation where h_c is the distance between flange centroids.
    h_c = h + tf
    
    e = (tf * h_c**2 * (b1**2 + b2**2)) / (8 * I_x_total)

    print("--- Shear Center Calculation ---")
    print("The shear center is located outside the cross-section, offset from the web.")
    print("Its distance 'e' from the web centerline is calculated based on the section's dimensions.")
    print("\nEquation (thin-wall approximation): e = (tf * h_c^2 * (b1^2 + b2^2)) / (8 * I_x)")
    
    # Print the equation with substituted values
    print("\nCalculation:")
    print(f"e = ({tf:.1f} * {h_c:.1f}^2 * ({b1:.1f}^2 + {b2:.1f}^2)) / (8 * {I_x_total:.2f})")
    
    # Print the intermediate steps
    numerator = tf * h_c**2 * (b1**2 + b2**2)
    denominator = 8 * I_x_total
    print(f"e = {numerator:.2f} / {denominator:.2f}")
    
    # Print the final result
    print(f"\nResult: The shear center is located e = {e:.2f} mm from the web centerline.")

calculate_shear_center_asymmetric_channel()