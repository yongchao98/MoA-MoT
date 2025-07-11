import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for a sample asymmetric channel section.
    The shear center is the point where a shear force can be applied without causing torsion.
    For a channel section, it is located at an offset 'e' from the web.
    """
    # --- Define Section Dimensions (in mm) ---
    # Asymmetric C-section with different top and bottom flanges
    B1 = 100.0  # Top flange total width
    tf1 = 15.0  # Top flange thickness
    B2 = 80.0   # Bottom flange total width
    tf2 = 12.0  # Bottom flange thickness
    H = 200.0   # Total height of the section
    tw = 10.0   # Web thickness

    print("--- Shear Center Calculation for an Asymmetric Channel Section ---\n")
    print(f"Dimensions (mm):")
    print(f"Top Flange: Width B1 = {B1}, Thickness tf1 = {tf1}")
    print(f"Bottom Flange: Width B2 = {B2}, Thickness tf2 = {tf2}")
    print(f"Total Height H = {H}, Web Thickness tw = {tw}\n")

    # --- 1. Calculate Centroid (y_bar) from the bottom edge ---
    h_web = H - tf1 - tf2 # Clear web height
    
    # Areas of components
    A_top = B1 * tf1
    A_web = h_web * tw
    A_bot = B2 * tf2
    A_total = A_top + A_web + A_bot

    # Y-coordinates of component centroids from bottom edge
    y_c_top = H - tf1 / 2.0
    y_c_web = tf2 + h_web / 2.0
    y_c_bot = tf2 / 2.0
    
    # Overall centroid (y_bar)
    y_bar = (A_top * y_c_top + A_web * y_c_web + A_bot * y_c_bot) / A_total
    
    print("--- Step 1: Locate the Centroid (Neutral Axis) ---")
    print(f"Area_top = {A_top:.2f} mm^2, Area_web = {A_web:.2f} mm^2, Area_bot = {A_bot:.2f} mm^2")
    print(f"y_bar = (A_top*y_c_top + A_web*y_c_web + A_bot*y_c_bot) / A_total")
    print(f"y_bar = ({A_top:.2f}*{y_c_top:.2f} + {A_web:.2f}*{y_c_web:.2f} + {A_bot:.2f}*{y_c_bot:.2f}) / {A_total:.2f}")
    print(f"Centroid location y_bar = {y_bar:.2f} mm (from the bottom edge)\n")

    # --- 2. Calculate Moment of Inertia (I_x) about the centroidal axis ---
    # Distances from component centroids to the overall centroid (for parallel axis theorem)
    d_top = y_c_top - y_bar
    d_web = y_c_web - y_bar
    d_bot = y_c_bot - y_bar
    
    # Local moments of inertia
    I_local_top = B1 * tf1**3 / 12.0
    I_local_web = tw * h_web**3 / 12.0
    I_local_bot = B2 * tf2**3 / 12.0
    
    # Total moment of inertia using Parallel Axis Theorem: I = I_local + A*d^2
    I_x = (I_local_top + A_top * d_top**2) + \
          (I_local_web + A_web * d_web**2) + \
          (I_local_bot + A_bot * d_bot**2)
          
    print("--- Step 2: Calculate Moment of Inertia (I_x) ---")
    print(f"I_x = (I_top_local + A_top*d_top^2) + (I_web_local + A_web*d_web^2) + (I_bot_local + A_bot*d_bot^2)")
    print(f"I_x = ({I_local_top:.2f} + {A_top:.2f}*{d_top:.2f}^2) + ({I_local_web:.2f} + {A_web:.2f}*{d_web:.2f}^2) + ({I_local_bot:.2f} + {A_bot:.2f}*{d_bot:.2f}^2)")
    print(f"Moment of Inertia I_x = {I_x:.2f} mm^4\n")

    # --- 3. Calculate Shear Center offset (e) from the web centerline ---
    # Using the formula for thin-walled open sections.
    # Flange widths measured from web face
    b1 = B1 - tw 
    b2 = B2 - tw

    # The integral part of the standard formula simplifies for this shape to:
    # integral(y*w*ds) where w is flange width from centerline, y is distance to NA
    # Note: Using a simplified but standard formula for thin-walled sections.
    term_top = (tf1 * d_top**2 * (B1**2 - tw**2)) / 4
    term_bot = (tf2 * d_bot**2 * (B2**2 - tw**2)) / 4
    
    # Using a more direct formula for offset 'e' from web centerline
    # e = (1 / I_x) * [Integral(omega*y*dA)] over flanges
    # For thin flanges, this approximates to:
    e_val = (1 / I_x) * ( (tf1*d_top*(B1**2/2 - tw**2/8)) + (tf2*abs(d_bot)*(B2**2/2 - tw**2/8)) )

    print("--- Step 3: Calculate Shear Center Offset (e) ---")
    print("Formula: e = (1/I_x) * [Integral(omega * y * dA)] for flanges")
    print("This evaluates to e = (1/I_x) * [tf1*d_top*(B1^2/2 - tw^2/8) + tf2*|d_bot|*(B2^2/2 - tw^2/8)]")
    
    # Numbers for the final equation
    p1 = tf1*d_top*(B1**2/2 - tw**2/8)
    p2 = tf2*abs(d_bot)*(B2**2/2 - tw**2/8)
    print(f"e = (1 / {I_x:.2f}) * [({tf1:.2f}*{d_top:.2f}*({B1**2/2:.2f} - {tw**2/8:.2f})) + ({tf2:.2f}*{abs(d_bot):.2f}*({B2**2/2:.2f} - {tw**2/8:.2f}))]")
    print(f"e = (1 / {I_x:.2f}) * [{p1:.2f} + {p2:.2f}]")
    print(f"Shear center offset e = {e_val:.2f} mm (from the web centerline)\n")
    
    print("Conclusion: The shear center is located outside the cross-section.")
    print("Its location is offset from the web by a distance 'e', which is determined by the channel's dimensions.")

if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()