import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for a thin-walled, asymmetric channel section.

    The shear center is the point through which shear loads must act to produce
    bending without twisting. For a channel section, it is located outside the
    section, on the side of the web.

    The plan is as follows:
    1. Define the geometry of the asymmetric channel section.
    2. Calculate the vertical location of the centroid (neutral axis), y_bar.
    3. Calculate the moment of inertia (I_x) about the horizontal neutral axis.
    4. Use the derived formula to calculate the shear center offset 'e' from the web.
    5. Print all calculations and the final result.
    """
    # 1. Define Geometry (all units in mm)
    b1 = 60   # Top flange width
    t1 = 10   # Top flange thickness
    b2 = 100  # Bottom flange width
    t2 = 10   # Bottom flange thickness
    h = 150   # Web height
    tw = 8    # Web thickness

    print("--- Section Properties ---")
    print(f"Top Flange: width(b1)={b1}, thickness(t1)={t1}")
    print(f"Bottom Flange: width(b2)={b2}, thickness(t2)={t2}")
    print(f"Web: height(h)={h}, thickness(tw)={tw}\n")

    # 2. Calculate Centroid (y_bar) from the bottom of the web
    A1 = b1 * t1  # Area of top flange
    y1 = h + t1 / 2.0  # y-centroid of top flange

    A_web = h * tw  # Area of web
    y_web = h / 2.0   # y-centroid of web

    A2 = b2 * t2  # Area of bottom flange
    y2 = -t2 / 2.0  # y-centroid of bottom flange

    total_area = A1 + A_web + A2
    sum_ay = (A1 * y1) + (A_web * y_web) + (A2 * y2)
    y_bar = sum_ay / total_area

    print("--- Centroid Calculation (from bottom of web) ---")
    print(f"y_bar = (A1*y1 + A_web*y_web + A2*y2) / (A1 + A_web + A2)")
    print(f"y_bar = ({A1:.0f}*{y1:.1f} + {A_web:.0f}*{y_web:.1f} + {A2:.0f}*{y2:.1f}) / ({A1:.0f} + {A_web:.0f} + {A2:.0f})")
    print(f"y_bar = {y_bar:.2f} mm\n")

    # 3. Calculate Moment of Inertia (I_x) using Parallel Axis Theorem
    # I_x = Σ (I_c + A * d^2)
    # d is the distance from the component's centroid to the section's centroid

    # Top Flange
    I_c1 = (b1 * t1**3) / 12.0
    d1 = y1 - y_bar
    I_x1 = I_c1 + A1 * d1**2

    # Web
    I_c_web = (tw * h**3) / 12.0
    d_web = y_web - y_bar
    I_x_web = I_c_web + A_web * d_web**2

    # Bottom Flange
    I_c2 = (b2 * t2**3) / 12.0
    d2 = y2 - y_bar
    I_x2 = I_c2 + A2 * d2**2

    I_x = I_x1 + I_x_web + I_x2

    print("--- Moment of Inertia Calculation (I_x) ---")
    print(f"I_x = I_x(top_flange) + I_x(web) + I_x(bottom_flange)")
    print(f"I_x = {I_x1:.2f} + {I_x_web:.2f} + {I_x2:.2f}")
    print(f"I_x = {I_x:.2f} mm^4\n")

    # 4. Calculate Shear Center offset 'e'
    # This formula calculates the moment produced by shear flow in the flanges.
    # e = (Moment from flange forces) / (Applied Shear Force V)
    # The V term cancels out, leaving a purely geometric formula.
    y_top_flange_center = h + t1/2 - y_bar # same as d1
    y_bottom_flange_center = y_bar - (-t2/2) # magnitude of d2

    # The formula is e = (1 / (2*I_x)) * (t1*b1^2*y_top^2 + t2*b2^2*y_bottom^2)
    # This formula contains a typo in many sources. A more robust formula derived
    # from balancing the torque from flange shear forces is used here.
    # Internal Torque T = (Force_top * h) + (Force_bottom * h), this is also wrong.
    # Let's use the correct formula:
    # Moment from top flange = (V / I_x) * (t1 * b1^2 / 2) * y_top_flange_center
    # Moment from bottom flange = (V / I_x) * (t2 * b2^2 / 2) * y_bottom_flange_center
    # The external moment V*e must balance the sum of the moments of the flange forces about the web.
    # The flange forces have a lever arm of h. Let's use a standard formula for simplicity
    # and correctness.
    e = (h**2 / (4 * I_x)) * (b1**2 * t1 + b2**2 * t2) # Note: this is an approximation for symmetric sections
    
    # Let's use a more accurate derivation for asymmetric section
    # e * V = Torque_top + Torque_bottom
    # Torque_top = Force_top * (h/2 + y_bar) ... no
    # Using a verified formula from a reliable source:
    # Let h_c = h + (t1+t2)/2 be distance between flange centroids
    h_c = h + (t1/2) + (t2/2)
    
    # Moment contribution from top flange (about NA)
    omega_1 = (h_c * b1 * t1) / 2
    # Moment contribution from bottom flange
    omega_2 = (h_c * b2 * t2) / 2
    
    # This is also getting complicated. Let's use the formula based on y distances, which is correct.
    y_top_center_from_NA = abs(y1 - y_bar)
    y_bot_center_from_NA = abs(y2 - y_bar)
    
    numerator = (t1 * b1**2 * y_top_center_from_NA**2) + (t2 * b2**2 * y_bot_center_from_NA**2)
    e = numerator / (2 * I_x)

    print("--- Shear Center Calculation (e) ---")
    print("The shear center offset 'e' from the web is calculated as:")
    print("e = (t1*b1²*y_top² + t2*b2²*y_bot²) / (2*I_x)")
    print(f"where y_top is distance from NA to top flange centroid = {y_top_center_from_NA:.2f}")
    print(f"and y_bottom is distance from NA to bottom flange centroid = {y_bot_center_from_NA:.2f}\n")
    
    print("Final Equation:")
    print(f"e = ({t1}*{b1}²*{y_top_center_from_NA:.2f}² + {t2}*{b2}²*{y_bot_center_from_NA:.2f}²) / (2*{I_x:.2f})")
    print(f"e = ({numerator:.2f}) / ({2*I_x:.2f})")
    print(f"e = {e:.2f} mm")
    print("\nThe shear center is located outside the section, at a distance 'e' from the web.")

calculate_shear_center_asymmetric_channel()