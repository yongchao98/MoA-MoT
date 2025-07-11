import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for a thin-walled, asymmetric channel section.
    The shear center is the point through which a shear force can be applied without causing torsion.
    For a channel section, it is located outside the section, at a distance 'e' from the web.
    """
    # --- 1. Define Section Geometry (all units in mm) ---
    # Top flange width
    b1 = 80.0
    # Bottom flange width
    b2 = 120.0
    # Web height (clear distance between flanges)
    hw = 200.0
    # Flange thickness
    tf = 10.0
    # Web thickness
    tw = 8.0

    print("--- Geometric Properties (mm) ---")
    print(f"Top Flange Width (b1):    {b1}")
    print(f"Bottom Flange Width (b2): {b2}")
    print(f"Web Height (hw):          {hw}")
    print(f"Flange Thickness (tf):    {tf}")
    print(f"Web Thickness (tw):       {tw}")
    print("-" * 35)

    # --- 2. Calculate Centroid (y_bar) Location ---
    # We set the origin (y=0) at the outside face of the bottom flange.
    # Area and y-centroid of each part relative to the origin:
    # Bottom Flange
    a_bot = b2 * tf
    y_bot = tf / 2.0
    # Web
    a_web = hw * tw
    y_web = tf + (hw / 2.0)
    # Top Flange
    a_top = b1 * tf
    y_top = tf + hw + (tf / 2.0)

    # Total Area
    total_area = a_top + a_web + a_bot
    # y_bar (Neutral Axis location from the bottom face)
    y_bar = (a_top * y_top + a_web * y_web + a_bot * y_bot) / total_area

    print(f"--- Centroid Calculation ---")
    print(f"Neutral Axis location (y_bar) from bottom face: {y_bar:.2f} mm")
    print("-" * 35)

    # --- 3. Calculate Moment of Inertia (I_x) about Centroidal Axis ---
    # Using the Parallel Axis Theorem: I_x = sum(I_local + A * d^2)
    # Distance from overall centroid (y_bar) to local centroid of each part
    d_top = y_top - y_bar
    d_web = y_web - y_bar
    d_bot = y_bot - y_bar

    # Local moments of inertia
    I_local_top = (b1 * tf**3) / 12.0
    I_local_web = (tw * hw**3) / 12.0
    I_local_bot = (b2 * tf**3) / 12.0

    # Total moment of inertia
    I_x = (I_local_top + a_top * d_top**2) + \
          (I_local_web + a_web * d_web**2) + \
          (I_local_bot + a_bot * d_bot**2)

    print(f"--- Moment of Inertia Calculation ---")
    print(f"Moment of Inertia (I_x): {I_x:,.2f} mm^4")
    print("-" * 35)

    # --- 4. Calculate Shear Center Eccentricity (e) ---
    # h1 and h2 are the distances from the neutral axis to the centerline of each flange
    h1 = abs(d_top)  # Distance to top flange centerline
    h2 = abs(d_bot)  # Distance to bottom flange centerline

    # Formula for shear center eccentricity from the web centerline:
    # e = (tf / (2 * I_x)) * (h1^2 * b1^2 + h2^2 * b2^2)
    e = (tf / (2 * I_x)) * (h1**2 * b1**2 + h2**2 * b2**2)
    
    print("--- Shear Center Calculation ---")
    print("The shear center eccentricity 'e' is calculated from the web's centerline.")
    print("Formula: e = (tf / (2 * I_x)) * (h1^2 * b1^2 + h2^2 * b2^2)\n")

    # Print the final equation with all numbers plugged in
    print("Final Equation:")
    print(f"e = ({tf} / (2 * {I_x:.2f})) * ({h1:.2f}\u00b2 * {b1}\u00b2 + {h2:.2f}\u00b2 * {b2}\u00b2)")

    # Print the final result
    print(f"\nResult:")
    print(f"The shear center is located at a distance e = {e:.2f} mm from the web.")
    print("-" * 35)

# Run the calculation
calculate_shear_center_asymmetric_channel()