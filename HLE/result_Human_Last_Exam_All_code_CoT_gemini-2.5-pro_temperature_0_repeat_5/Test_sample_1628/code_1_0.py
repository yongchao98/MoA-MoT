import math

def calculate_shear_center():
    """
    Calculates the location of the shear center for an asymmetric channel section.

    The shear center is a point where a shear force can be applied without causing torsion.
    For an asymmetric channel, this point is located outside the section, offset from the web.
    This script demonstrates the calculation of this offset based on the section's dimensions.
    """
    # --- Step 1: Define the geometry of the asymmetric channel section ---
    # Using centerline dimensions in millimeters (mm).
    h = 150.0     # Height of the web (distance between flange centerlines)
    t_w = 6.0     # Thickness of the web
    b1 = 100.0    # Width of the top flange (measured from web centerline)
    t_f1 = 10.0   # Thickness of the top flange
    b2 = 50.0     # Width of the bottom flange (measured from web centerline)
    t_f2 = 8.0    # Thickness of the bottom flange

    print("--- Input Section Properties ---")
    print(f"Web height (h): {h} mm")
    print(f"Web thickness (t_w): {t_w} mm")
    print(f"Top flange width (b1): {b1} mm, thickness (t_f1): {t_f1} mm")
    print(f"Bottom flange width (b2): {b2} mm, thickness (t_f2): {t_f2} mm")

    # --- Step 2: Calculate the vertical location of the centroid (y_bar) ---
    # Origin is at the web's mid-height.
    A_w = h * t_w
    A1 = b1 * t_f1
    A2 = b2 * t_f2
    A_total = A_w + A1 + A2
    
    # y_bar = (Sum of A_i * y_i) / A_total
    y_bar = (A1 * (h / 2) + A2 * (-h / 2)) / A_total

    # --- Step 3: Calculate the Moment of Inertia (I_x) about the centroidal x-axis ---
    # Using the Parallel Axis Theorem: I = I_c + A*d^2
    d_w_sq = (-y_bar)**2
    d1_sq = (h / 2 - y_bar)**2
    d2_sq = (-h / 2 - y_bar)**2

    I_xw = (t_w * h**3) / 12 + A_w * d_w_sq
    I_x1 = (b1 * t_f1**3) / 12 + A1 * d1_sq
    I_x2 = (b2 * t_f2**3) / 12 + A2 * d2_sq
    I_x = I_xw + I_x1 + I_x2

    # --- Step 4: Calculate the Shear Center horizontal eccentricity (e) ---
    # The formula for the horizontal offset 'e' from the web's centerline is:
    # e = (t_f1*b1^2*d1^2 + t_f2*b2^2*d2^2) / (2*I_x)
    # where d1 and d2 are the squared distances from the neutral axis to the flange centerlines.
    
    numerator = (t_f1 * b1**2 * d1_sq + t_f2 * b2**2 * d2_sq)
    denominator = 2 * I_x
    e = numerator / denominator

    print("\n--- Calculation of Shear Center Location ---")
    print("The shear center's horizontal offset 'e' from the web is calculated as follows:")
    print("\n1. Vertical centroid from web mid-height (y_bar): {:.2f} mm".format(y_bar))
    print("2. Moment of Inertia about neutral axis (I_x): {:.2f} mm^4".format(I_x))
    print("3. Distances from neutral axis to flange centroids squared:")
    print("   d1^2 = {:.2f} mm^2, d2^2 = {:.2f} mm^2".format(d1_sq, d2_sq))

    print("\nFinal Equation for eccentricity 'e':")
    print("e = (t_f1*b1^2*d1^2 + t_f2*b2^2*d2^2) / (2*I_x)")
    print("e = ({:.1f}*{:.1f}^2*{:.2f} + {:.1f}*{:.1f}^2*{:.2f}) / (2*{:.2f})".format(
        t_f1, b1, d1_sq, t_f2, b2, d2_sq, I_x))
    print("e = {:.2f} / {:.2f}".format(numerator, denominator))
    print("\nResult:")
    print("The shear center is located {:.2f} mm outside the web's centerline.".format(e))

if __name__ == '__main__':
    calculate_shear_center()