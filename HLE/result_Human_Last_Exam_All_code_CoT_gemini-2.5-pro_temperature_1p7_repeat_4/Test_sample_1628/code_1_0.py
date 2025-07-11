import math

def calculate_asymmetric_channel_shear_center(b1, t1, b2, t2, h_total, tw):
    """
    Calculates the shear center for an asymmetric channel section.
    
    The shear center 'e' is the horizontal distance from the web's centerline.
    
    Args:
        b1 (float): Width of the top flange.
        t1 (float): Thickness of the top flange.
        b2 (float): Width of the bottom flange.
        t2 (float): Thickness of the bottom flange.
        h_total (float): Total height of the section.
        tw (float): Thickness of the web.
    """
    
    # 1. Define component areas and their local centroids (from the bottom edge of the section)
    # Flange 1 (top)
    a1 = b1 * t1
    y1 = h_total - t1 / 2.0
    
    # Flange 2 (bottom)
    a2 = b2 * t2
    y2 = t2 / 2.0
    
    # Web
    # Effective web height between the inside faces of the flanges
    h_web = h_total - t1 - t2
    aw = h_web * tw
    yw = t2 + h_web / 2.0
    
    # Total area
    A_total = a1 + a2 + aw
    
    # 2. Calculate the vertical centroid (y_bar) from the bottom edge
    y_bar = (a1 * y1 + a2 * y2 + aw * yw) / A_total
    
    # 3. Calculate the moment of inertia (Ix) about the centroidal x-axis (parallel axis theorem)
    # Ix = I_local + A*d^2
    # I_local for a rectangle = (base * height^3) / 12
    
    # Ix for Flange 1 (top)
    d1 = y1 - y_bar
    Ix1 = (b1 * t1**3) / 12.0 + a1 * d1**2
    
    # Ix for Flange 2 (bottom)
    d2 = y2 - y_bar
    Ix2 = (b2 * t2**3) / 12.0 + a2 * d2**2
    
    # Ix for Web
    dw = yw - y_bar
    Ix_w = (tw * h_web**3) / 12.0 + aw * dw**2
    
    Ix_total = Ix1 + Ix2 + Ix_w
    
    # 4. Calculate shear center eccentricity 'e' from the web centerline
    # h_c1 and h_c2 are the distances from the neutral axis to the flange centerlines
    h_c1 = y1 - y_bar
    h_c2 = abs(y2 - y_bar) # Use absolute value
    
    # b1_out and b2_out are flange lengths from the web face
    b1_out = b1 - tw/2.0
    b2_out = b2 - tw/2.0
    
    # Formula for shear center eccentricity 'e'
    # Derived from balancing the internal moment from flange shear flows
    term1 = (t1 * (h_c1**2) * (b1_out**2)) / 2.0
    term2 = (t2 * (h_c2**2) * (b2_out**2)) / 2.0
    
    e = (term1 + term2) / Ix_total
    
    # 5. Print results
    print("--- Asymmetric Channel Section Properties ---")
    print(f"Top Flange (b1 x t1):      {b1} x {t1}")
    print(f"Bottom Flange (b2 x t2):   {b2} x {t2}")
    print(f"Web (h_total x tw):      {h_total} x {tw}\n")
    
    print("--- Calculated Intermediate Values ---")
    print(f"Centroid location (y_bar from bottom): {y_bar:.2f}")
    print(f"Moment of Inertia (Ix):                {Ix_total:.2f}\n")

    print("--- Final Result ---")
    print("The shear center is a point through which applied shear forces produce no twisting.")
    print("For this channel section, it is located outside the section.")
    print(f"Calculated distance 'e' from web centerline: {e:.2f}")
    print("\nThis demonstrates that the shear center is located outside the cross-section, offset from the centroid, at a distance determined by the channel's dimensions.")

# Example Dimensions (e.g., in mm)
top_flange_width = 100
top_flange_thickness = 10
bottom_flange_width = 150
bottom_flange_thickness = 12
total_height = 200
web_thickness = 8

calculate_asymmetric_channel_shear_center(top_flange_width, top_flange_thickness,
                                          bottom_flange_width, bottom_flange_thickness,
                                          total_height, web_thickness)