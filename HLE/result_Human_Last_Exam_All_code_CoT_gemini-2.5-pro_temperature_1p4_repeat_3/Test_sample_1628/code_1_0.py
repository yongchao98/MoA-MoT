import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for an asymmetric channel section.

    The shear center is the point through which an applied shear force
    produces no twisting moment. For an asymmetric channel, this point
    is located outside the section, at an eccentricity 'e' from the web.
    """
    # 1. Define Geometry (in mm)
    h = 150.0   # Web height
    b1 = 80.0   # Top flange width
    b2 = 60.0   # Bottom flange width
    t_w = 6.0   # Web thickness
    t_f1 = 8.0  # Top flange thickness
    t_f2 = 10.0 # Bottom flange thickness

    print("--- Section Dimensions ---")
    print(f"Web Height (h): {h} mm")
    print(f"Top Flange Width (b1): {b1} mm")
    print(f"Bottom Flange Width (b2): {b2} mm")
    print(f"Web Thickness (t_w): {t_w} mm")
    print(f"Top Flange Thickness (t_f1): {t_f1} mm")
    print(f"Bottom Flange Thickness (t_f2): {t_f2} mm")
    print("-" * 28 + "\n")

    # 2. Calculate Centroid (y_c)
    # The origin (y=0) is the bottom face of the bottom flange.
    # Areas of components
    A1 = b1 * t_f1  # Top flange
    A2 = h * t_w    # Web
    A3 = b2 * t_f2  # Bottom flange
    A_total = A1 + A2 + A3

    # y-coordinates of component centroids from the origin
    y1_comp = t_f2 + h - t_f1 / 2.0  # Centroid of top flange
    y2_comp = t_f2 + h / 2.0         # Centroid of web
    y3_comp = t_f2 / 2.0           # Centroid of bottom flange

    # Vertical centroid (y_c)
    y_c = (A1 * y1_comp + A2 * y2_comp + A3 * y3_comp) / A_total

    print("--- Centroid and Moment of Inertia Calculation ---")
    print(f"Vertical Centroid (y_c from bottom): {y_c:.2f} mm")

    # 3. Calculate Moment of Inertia (Ix) about the Centroidal Axis
    # Using the Parallel Axis Theorem: Ix = Σ(I_local + A * d^2)
    I1_local = b1 * t_f1**3 / 12.0
    I2_local = t_w * h**3 / 12.0
    I3_local = b2 * t_f2**3 / 12.0
    
    d1 = y1_comp - y_c
    d2 = y2_comp - y_c
    d3 = y3_comp - y_c

    Ix = (I1_local + A1 * d1**2) + (I2_local + A2 * d2**2) + (I3_local + A3 * d3**2)
    print(f"Moment of Inertia (Ix): {Ix:.2f} mm^4")
    print("-" * 50 + "\n")

    # 4. Calculate Shear Center Eccentricity (e)
    # 'e' is the offset of the shear center from the centerline of the web.
    # h1 and h2 are the distances from the centroid to the centerlines of the flanges.
    h1 = abs(d1)
    h2 = abs(d3)
    
    print("--- Shear Center Calculation ---")
    print(f"Distance from centroid to top flange centerline (h1): {h1:.2f} mm")
    print(f"Distance from centroid to bottom flange centerline (h2): {h2:.2f} mm\n")

    # The formula for eccentricity 'e' is:
    # e = (t_f1*b1^2*h1^2 + t_f2*b2^2*h2^2) / (2 * Ix)
    numerator = (t_f1 * b1**2 * h1**2) + (t_f2 * b2**2 * h2**2)
    denominator = 2 * Ix
    e = numerator / denominator

    # Output the final equation with numerical values
    print("Formula for shear center eccentricity 'e':")
    print("e = (t_f1 * b1^2 * h1^2 + t_f2 * b2^2 * h2^2) / (2 * Ix)\n")
    print("Substituting the values into the formula:")
    print(f"e = ({t_f1:.2f}*{b1:.2f}^2*{h1:.2f}^2 + {t_f2:.2f}*{b2:.2f}^2*{h2:.2f}^2) / (2 * {Ix:.2f})")
    print(f"e = ({t_f1 * b1**2 * h1**2:.2f} + {t_f2 * b2**2 * h2**2:.2f}) / ({denominator:.2f})")
    print(f"e = {numerator:.2f} / {denominator:.2f}")
    
    print("\n--- Final Result ---")
    print(f"The shear center is located at an eccentricity (e) of {e:.2f} mm from the web's centerline.")
    print("This point is outside the cross-section, ensuring that a shear force applied here does not cause twisting.")

if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()