import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for a defined asymmetric channel section.

    The shear center is the point through which shear loads must act to produce
    no twisting. For a channel section, it is located outside the section,
    offset from the web.

    This script calculates the location by:
    1. Defining the section's geometry.
    2. Calculating the area and centroid of the section.
    3. Calculating the moment of inertia (I_x) about the centroidal x-axis.
    4. Using a formula to find the horizontal offset 'e' of the shear center
       from the centroidal vertical axis.
    """
    # 1. Define Section Geometry (all units in mm)
    b1 = 80.0   # Top flange width
    t1 = 10.0   # Top flange thickness
    b2 = 50.0   # Bottom flange width
    t2 = 10.0   # Bottom flange thickness
    h_web = 130.0 # Height of the web element
    tw = 8.0    # Web thickness
    h_total = h_web + t1 + t2

    print("--- Section Properties (mm) ---")
    print(f"Top Flange (b1 x t1):     {b1} x {t1}")
    print(f"Bottom Flange (b2 x t2):    {b2} x {t2}")
    print(f"Web (h_web x tw):         {h_web} x {tw}")
    print(f"Total Height:               {h_total}")
    print("-" * 33)

    # 2. Calculate Area and Centroid (y_c from the bottom edge)
    # The section is divided into three rectangular parts:
    # 1: Bottom flange, 2: Web, 3: Top flange
    A1 = b2 * t2
    y1 = t2 / 2
    A2 = h_web * tw
    y2 = t2 + h_web / 2
    A3 = b1 * t1
    y3 = t2 + h_web + t1 / 2

    A_total = A1 + A2 + A3
    y_c = (A1 * y1 + A2 * y2 + A3 * y3) / A_total

    print("--- Calculated Properties ---")
    print(f"Total Area:                 {A_total:.2f} mm^2")
    print(f"Centroid (y_c from bottom): {y_c:.2f} mm")

    # 3. Calculate Moment of Inertia (I_x) about the centroidal x-axis
    # Using the Parallel Axis Theorem: I_x = sum(I_ci + A_i * d_i^2)
    # I_ci for a rectangle is (b*h^3)/12
    I_c1 = (b2 * t2**3) / 12
    d_c1 = y_c - y1
    I_c2 = (tw * h_web**3) / 12
    d_c2 = y_c - y2
    I_c3 = (b1 * t1**3) / 12
    d_c3 = y3 - y_c
    
    I_x = (I_c1 + A1 * d_c1**2) + (I_c2 + A2 * d_c2**2) + (I_c3 + A3 * d_c3**2)
    print(f"Moment of Inertia (I_x):    {I_x:.2f} mm^4")

    # 4. Calculate Shear Center offset 'e' from the centroidal vertical axis
    # The offset 'e' is given by the formula which balances the external shear
    # force moment with the internal twisting moment from the shear flow in flanges.
    # Formula: e = (1 / (2*I_x)) * (t1*d1^2*b1^2 + t2*d2^2*b2^2)
    # where d1 and d2 are the distances from the section centroid to the
    # centroids of the top and bottom flanges, respectively.
    d1_flange = d_c3 # This is the distance to the top flange centroid
    d2_flange = d_c1 # This is the distance to the bottom flange centroid
    
    numerator = (t1 * d1_flange**2 * b1**2) + (t2 * d2_flange**2 * b2**2)
    e = numerator / (2 * I_x)

    print("\n--- Final Result ---")
    # Note: For a typical channel shape, the centroid's x-coordinate is inside
    # the flanges, and this calculated offset 'e' will place the shear center
    # outside the material section, to the left of the web.
    print(f"Horizontal offset of the Shear Center from the Centroid (e): {e:.2f} mm")
    print("\nThis result demonstrates that the shear center is located outside the cross-section,")
    print("offset from the web, at a distance calculated based on the section's dimensions.")

if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()