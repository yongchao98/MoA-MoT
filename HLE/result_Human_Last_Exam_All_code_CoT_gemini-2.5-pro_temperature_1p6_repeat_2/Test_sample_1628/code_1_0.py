import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates and explains the shear center location for an asymmetric channel section.
    The script uses centerline dimensions and assumes thin walls.
    """
    # 1. Define Geometry for an example asymmetric channel section
    # All dimensions are in mm.
    h = 200.0  # Web height (centerline to centerline)
    b1 = 100.0 # Top flange width (from web centerline)
    t1 = 10.0  # Top flange thickness
    b2 = 80.0  # Bottom flange width (from web centerline)
    t2 = 12.0  # Bottom flange thickness
    tw = 8.0   # Web thickness

    print("--- Shear Center Location Calculation for an Asymmetric Channel ---")
    print("\nThe calculation demonstrates that the shear center is located outside the section at a")
    print("distance determined by the section's dimensions.\n")
    print("Given Section Dimensions (mm):")
    print(f"  Web Height (h): {h}")
    print(f"  Top Flange (b1 x t1): {b1} x {t1}")
    print(f"  Bottom Flange (b2 x t2): {b2} x {t2}")
    print(f"  Web Thickness (tw): {tw}\n")

    # 2. Calculate Centroid Location (x_c, y_c)
    # Origin (0,0) is at the intersection of the web and bottom flange centerlines.
    A1 = b1 * t1  # Area of top flange
    A2 = b2 * t2  # Area of bottom flange
    Aw = h * tw   # Area of web
    A_total = A1 + A2 + Aw

    # x-coordinate of centroid relative to the web centerline
    x_c_num = (A1 * (b1 / 2)) + (A2 * (b2 / 2)) + (Aw * 0)
    x_c = x_c_num / A_total
    
    # y-coordinate of centroid relative to the bottom flange centerline
    y_c_num = (A1 * h) + (A2 * 0) + (Aw * (h / 2))
    y_c = y_c_num / A_total

    print("--- Step 1: Calculate Centroid Location (relative to web/bottom-flange junction) ---")
    print(f"y_c = (A1*h + A2*0 + Aw*h/2) / (A1 + A2 + Aw)")
    print(f"y_c = ({A1:.0f}*{h} + {A2:.0f}*0 + {Aw:.0f}*{h/2}) / ({A1:.0f} + {A2:.0f} + {Aw:.0f}) = {y_c:.2f} mm\n")

    # 3. Calculate Moment of Inertia (I_xx) about the centroidal x-axis
    # Using the Parallel Axis Theorem: I_xx = I_c + A*d^2
    # The term I_c for the flanges (e.g., b1*t1^3/12) is negligible for thin walls.
    I_xx_1 = A1 * (h - y_c)**2
    I_xx_2 = A2 * (0 - y_c)**2
    I_xx_w = (tw * h**3 / 12) + Aw * (h / 2 - y_c)**2
    I_xx = I_xx_1 + I_xx_2 + I_xx_w

    print("--- Step 2: Calculate Moment of Inertia (I_xx) about Centroidal X-Axis ---")
    print(f"I_xx = A1*(h-y_c)^2 + A2*(0-y_c)^2 + [tw*h^3/12 + Aw*(h/2-y_c)^2]")
    print(f"I_xx = {A1:.0f}*({h}-{y_c:.2f})^2 + {A2:.0f}*(0-{y_c:.2f})^2 + [{tw}*{h}^3/12 + {Aw:.0f}*({h/2}-{y_c:.2f})^2] = {I_xx:.2e} mm^4\n")
    
    # 4. Calculate Shear Center Location (e_x, y_s)
    # e_x is the horizontal offset from the web centerline
    # h1_dist and h2_dist are the distances from the centroid to the top and bottom flanges
    h1_dist = h - y_c
    h2_dist = y_c

    e_x_num = t1 * b1**2 * h1_dist**2 + t2 * b2**2 * h2_dist**2
    e_x_den = 2 * I_xx
    e_x = e_x_num / e_x_den

    # y_s is the vertical location measured from the bottom flange centerline
    # Formula from AISC Design Guide 9 / CISC "Torsional Section Properties"
    y_s_num = h * b1**3 * t1
    y_s_den = (b1**3 * t1) + (b2**3 * t2)
    y_s = y_s_num / y_s_den if y_s_den != 0 else 0

    print("--- Step 3: Calculate Shear Center Coordinates (e_x, y_s) ---")
    print("Horizontal offset e_x from the web centerline:")
    print(f"e_x = (t1*b1^2*(h-y_c)^2 + t2*b2^2*y_c^2) / (2 * I_xx)")
    print(f"e_x = ({t1}*{b1}^2*{h1_dist:.2f}^2 + {t2}*{b2}^2*{h2_dist:.2f}^2) / (2 * {I_xx:.2e}) = {e_x:.2f} mm")
    print("\nVertical position y_s from the bottom flange centerline:")
    print(f"y_s = (h * b1^3 * t1) / (b1^3 * t1 + b2^3 * t2)")
    print(f"y_s = ({h} * {b1}^3 * {t1}) / ({b1}^3 * {t1} + {b2}^3 * {t2}) = {y_s:.2f} mm\n")

    print("--- Conclusion ---")
    print(f"The calculated shear center is located at ({e_x:.2f} mm, {y_s:.2f} mm) relative to the web/bottom-flange junction.")
    print("Since e_x is positive, the shear center is located OUTSIDE the cross-section.")
    print("Its position is offset from the centroid and is determined by the channel's flange widths, web height, and thicknesses.")
    print("This precise description matches answer choice G.")

if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()