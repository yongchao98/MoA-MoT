import math

def calculate_shear_center_asymmetric_channel(b1, t1, b2, t2, h, tw):
    """
    Calculates the shear center for an asymmetric channel section.

    The origin (0,0) is assumed to be at the bottom-left corner of the web.
    The web lies along the y-axis from y=0 to y=h.
    The top flange extends from the web in the positive x-direction.
    The bottom flange extends from the web in the positive x-direction.

    Args:
        b1 (float): Width of the top flange.
        t1 (float): Thickness of the top flange.
        b2 (float): Width of the bottom flange.
        t2 (float): Thickness of the bottom flange.
        h (float): Clear height of the web.
        tw (float): Thickness of the web.
    """

    print("--- Input Section Properties ---")
    print(f"Top Flange (b1 x t1):       {b1:.2f} x {t1:.2f}")
    print(f"Bottom Flange (b2 x t2):    {b2:.2f} x {t2:.2f}")
    print(f"Web (h x tw):               {h:.2f} x {tw:.2f}\n")

    # 1. Calculate properties of each component (Area, Centroid coordinates)
    # Component 1: Top Flange
    a1 = b1 * t1
    x1 = b1 / 2
    y1 = h + t1 / 2

    # Component 2: Web
    a2 = h * tw
    x2 = -tw / 2  # Assuming web centerline is on y-axis for symmetry in x-calculation later
    y2 = h / 2
    
    # Component 3: Bottom Flange
    a3 = b2 * t2
    x3 = b2 / 2
    y3 = -t2 / 2

    # 2. Calculate the centroid of the entire section (x_bar, y_bar)
    # Origin for this calculation is the centerline of the web and the bottom edge of the web.
    total_area = a1 + a2 + a3
    # Y-centroid (y_bar) from bottom of web
    y_bar = (a1 * (h + t1/2) + a2 * (h/2) + a3 * (-t2/2)) / total_area
    # X-centroid (x_bar) from web centerline
    x_bar = (a1 * (tw/2 + b1/2) + a2 * 0 + a3 * (tw/2 + b2/2)) / total_area

    # 3. Calculate Moment of Inertia (Ix) about the horizontal centroidal axis
    # Using Parallel Axis Theorem: I = I_local + A*d^2
    # Distance from component centroid to overall y_bar
    d1_y = (h + t1/2) - y_bar
    d2_y = (h/2) - y_bar
    d3_y = (-t2/2) - y_bar
    
    i_c1 = (b1 * t1**3) / 12
    i_c2 = (tw * h**3) / 12
    i_c3 = (b2 * t2**3) / 12
    
    ix = (i_c1 + a1 * d1_y**2) + (i_c2 + a2 * d2_y**2) + (i_c3 + a3 * d3_y**2)
    
    # 4. Calculate the shear center location 'e'
    # The formula for the distance 'e0' of the shear center from the web's centerline is:
    # e0 = (1/Ix) * Integral(omega * y * dA), where omega is the sectorial coordinate.
    # For a channel section, this simplifies to:
    # h1 is distance from NA to top flange centerline, h2 is distance from NA to bottom flange centerline.
    # Note: These are absolute distances.
    h1 = abs(d1_y)
    h2 = abs(d3_y)

    # Note: For this formula, b1 and b2 are the flange widths extending from the web.
    e_offset_from_web_centerline = (t1 * b1**2 * h1 / 2 + t2 * b2**2 * h2 / 2) / ix
    
    print("--- Calculated Section Properties ---")
    print(f"Total Area:                 {total_area:.2f}")
    print(f"Y-Centroid (from web bottom): {y_bar:.2f}")
    print(f"Moment of Inertia (Ix):     {ix:.2f}\n")
    
    print("--- Shear Center Location ---")
    print("The shear center is a point offset horizontally from the web.")
    print(f"Calculation: e = (t1*b1^2*h1/2 + t2*b2^2*h2/2) / Ix")
    print(f"             e = ({t1:.2f}*{b1:.2f}^2*{h1:.2f}/2 + {t2:.2f}*{b2:.2f}^2*{h2:.2f}/2) / {ix:.2f}")
    print(f"Result: e = {e_offset_from_web_centerline:.2f}")
    print(f"\nThe shear center is located {e_offset_from_web_centerline:.2f} units from the centerline of the web, outside the section.")


# --- Example Asymmetric Channel Section ---
# All units are consistent (e.g., mm)
top_flange_width = 80
top_flange_thickness = 10
bottom_flange_width = 60
bottom_flange_thickness = 8
web_height = 150
web_thickness = 5

calculate_shear_center_asymmetric_channel(
    b1=top_flange_width, 
    t1=top_flange_thickness,
    b2=bottom_flange_width,
    t2=bottom_flange_thickness,
    h=web_height,
    tw=web_thickness
)