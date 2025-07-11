import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center for a thin-walled, asymmetric channel section.
    The shear center is the point through which an applied shear force
    produces no twisting moment. For a channel, this function calculates
    the horizontal offset 'e' from the centerline of the web.
    """
    # --- 1. Define Section Dimensions (in consistent units, e.g., mm) ---
    # Example dimensions for an asymmetric channel section
    b1 = 150.0  # Top flange total width
    tf1 = 12.0  # Top flange thickness
    b2 = 100.0  # Bottom flange total width
    tf2 = 10.0  # Bottom flange thickness
    h_w = 200.0 # Web height (between flanges)
    t_w = 8.0   # Web thickness

    print("--- Input Section Dimensions (mm) ---")
    print(f"Top Flange (b1 x tf1):      {b1} x {tf1}")
    print(f"Bottom Flange (b2 x tf2):   {b2} x {tf2}")
    print(f"Web (h_w x t_w):            {h_w} x {t_w}\n")

    # --- 2. Calculate Centroid (Neutral Axis y_bar) ---
    # We measure y from the outer face of the bottom flange.
    
    # Areas of components
    A1 = b1 * tf1  # Top flange
    A2 = h_w * t_w # Web
    A3 = b2 * tf2  # Bottom flange
    A_total = A1 + A2 + A3

    # Centroid of each component (y_i)
    y1_c = tf2 + h_w + tf1 / 2.0  # Centroid of top flange
    y2_c = tf2 + h_w / 2.0         # Centroid of web
    y3_c = tf2 / 2.0               # Centroid of bottom flange

    # Overall centroid (y_bar)
    y_bar = (A1 * y1_c + A2 * y2_c + A3 * y3_c) / A_total

    print("--- Step 1: Centroid Calculation ---")
    print(f"The vertical position of the centroid (y_bar) from the bottom edge is: {y_bar:.2f} mm\n")

    # --- 3. Calculate Moment of Inertia (I_x) about the centroidal x-axis ---
    # Using the Parallel Axis Theorem: I = I_local + A*d^2
    
    # Moment of inertia for each component about its own centroid
    I1_local = (b1 * tf1**3) / 12.0
    I2_local = (t_w * h_w**3) / 12.0
    I3_local = (b2 * tf2**3) / 12.0
    
    # Distances from overall centroid to component centroids
    d1 = y1_c - y_bar
    d2 = y2_c - y_bar
    d3 = y3_c - y_bar

    # Total moment of inertia using Parallel Axis Theorem
    I_x = (I1_local + A1 * d1**2) + \
          (I2_local + A2 * d2**2) + \
          (I3_local + A3 * d3**2)
          
    print("--- Step 2: Moment of Inertia Calculation ---")
    print(f"The moment of inertia (I_x) about the horizontal centroidal axis is: {I_x:,.2f} mm^4\n")
    
    # --- 4. Calculate Shear Center Offset (e) from the web centerline ---
    # The formula is derived from balancing the internal twisting moment from
    # flange shear flows with the external moment V*e.
    # e = (1 / (2*I_x)) * [t_f1 * b_o1^2 * d_1^2 + t_f2 * b_o2^2 * d_2^2]
    # where b_o is the flange outstand from the web face.
    
    b_o1 = b1 - t_w / 2.0  # Top flange outstand from web centerline
    b_o2 = b2 - t_w / 2.0  # Bottom flange outstand from web centerline

    # The distances d1 and d3 are the distances from the neutral axis to the
    # centerlines of the top and bottom flanges, respectively.
    d_top_flange = abs(d1)
    d_bot_flange = abs(d3)
    
    numerator = (tf1 * (b_o1**2) * (d_top_flange**2)) + (tf2 * (b_o2**2) * (d_bot_flange**2))
    e = numerator / (2 * I_x)
    
    print("--- Step 3: Shear Center Calculation ---")
    print(f"Formula for shear center offset 'e' from web centerline:")
    print("e = (tf1*b_o1^2*d_top^2 + tf2*b_o2^2*d_bot^2) / (2*I_x)")
    print("where: ")
    print(f"  tf1 = {tf1}, b_o1 (top flange outstand) = {b_o1}, d_top (top flange dist to NA) = {d_top_flange:.2f}")
    print(f"  tf2 = {tf2}, b_o2 (bot flange outstand) = {b_o2}, d_bot (bot flange dist to NA) = {d_bot_flange:.2f}")
    print(f"  I_x = {I_x:,.2f}\n")
    print("--- Final Result ---")
    print(f"The shear center is located at a horizontal distance 'e' = {e:.2f} mm from the centerline of the web.")
    print("This point is outside the physical cross-section.")

if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()