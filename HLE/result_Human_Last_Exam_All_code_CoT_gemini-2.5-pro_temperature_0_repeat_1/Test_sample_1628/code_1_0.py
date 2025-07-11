import math

def calculate_shear_center():
    """
    Calculates the shear center for a thin-walled, asymmetric channel section.
    """
    # --- 1. Define Section Geometry (in mm) ---
    # b1: Top flange width (from web centerline)
    # b2: Bottom flange width (from web centerline)
    # h: Web height (between flange centerlines)
    # t_f: Flange thickness
    # t_w: Web thickness
    b1 = 100.0
    b2 = 50.0
    h = 200.0
    t_f = 10.0
    t_w = 8.0

    print("--- Asymmetric Channel Section Properties ---")
    print(f"Top flange width (b1): {b1} mm")
    print(f"Bottom flange width (b2): {b2} mm")
    print(f"Web height (h): {h} mm")
    print(f"Flange thickness (t_f): {t_f} mm")
    print(f"Web thickness (t_w): {t_w} mm\n")

    # --- 2. Calculate Centroid (y_bar) ---
    # Using thin-walled approximation. Origin is at the center of the web.
    # y-coordinates of flange centerlines from web center:
    y_f1 = h / 2
    y_f2 = -h / 2

    # Areas of the components
    A_f1 = b1 * t_f  # Area of top flange
    A_f2 = b2 * t_f  # Area of bottom flange
    A_w = h * t_w    # Area of web
    A_total = A_f1 + A_f2 + A_w

    # Vertical centroid location (y_bar) from the web's horizontal centerline
    y_bar = (A_f1 * y_f1 + A_f2 * y_f2 + A_w * 0) / A_total

    print("--- Centroid Calculation (y_bar) ---")
    print(f"y_bar = (A_f1*y_f1 + A_f2*y_f2 + A_w*0) / (A_f1 + A_f2 + A_w)")
    print(f"y_bar = ({A_f1:.1f}*{y_f1:.1f} + {A_f2:.1f}*{y_f2:.1f} + {A_w:.1f}*0) / ({A_f1:.1f} + {A_f2:.1f} + {A_w:.1f})")
    print(f"y_bar = {y_bar:.2f} mm (from web's horizontal centerline)\n")

    # --- 3. Calculate Moment of Inertia (I_x) ---
    # Using Parallel Axis Theorem: I = I_local + A*d^2
    # I_local for thin flanges about their own axis is negligible.
    I_f1 = A_f1 * (y_f1 - y_bar)**2
    I_f2 = A_f2 * (y_f2 - y_bar)**2
    I_w = (t_w * h**3) / 12 + A_w * (0 - y_bar)**2
    I_x = I_f1 + I_f2 + I_w

    print("--- Moment of Inertia Calculation (I_x) ---")
    print("I_x = I_top_flange + I_bottom_flange + I_web")
    print(f"I_top_flange = A_f1 * (y_f1 - y_bar)^2 = {A_f1:.1f} * ({y_f1:.1f} - {y_bar:.2f})^2 = {I_f1:.2f} mm^4")
    print(f"I_bottom_flange = A_f2 * (y_f2 - y_bar)^2 = {A_f2:.1f} * ({y_f2:.1f} - {y_bar:.2f})^2 = {I_f2:.2f} mm^4")
    print(f"I_web = (t_w*h^3)/12 + A_w*(0 - y_bar)^2 = ({t_w:.1f}*{h:.1f}^3)/12 + {A_w:.1f}*(0 - {y_bar:.2f})^2 = {I_w:.2f} mm^4")
    print(f"I_x = {I_f1:.2f} + {I_f2:.2f} + {I_w:.2f} = {I_x:.2f} mm^4\n")

    # --- 4. Calculate Shear Center Offset (e) ---
    # e is the horizontal offset from the web's vertical centerline.
    # This formula sums the moments produced by the shear forces in the flanges.
    term1 = b1**2 * (y_f1 - y_bar)**2
    term2 = b2**2 * (y_f2 - y_bar)**2
    e = (t_f / (2 * I_x)) * (term1 + term2)

    print("--- Shear Center Offset Calculation (e) ---")
    print("e = (t_f / (2 * I_x)) * [b1^2 * (y_f1 - y_bar)^2 + b2^2 * (y_f2 - y_bar)^2]")
    print(f"e = ({t_f:.1f} / (2 * {I_x:.2f})) * [{b1**2:.1f} * ({y_f1:.1f} - {y_bar:.2f})^2 + {b2**2:.1f} * ({y_f2:.1f} - {y_bar:.2f})^2]")
    print(f"e = ({t_f / (2 * I_x):.4e}) * [{term1:.2f} + {term2:.2f}]")
    print(f"Final Result: e = {e:.2f} mm\n")

    print("--- Conclusion ---")
    print(f"The shear center is located {e:.2f} mm from the web's centerline, outside the section.")
    print("This confirms that for an asymmetric channel, the shear center is offset from the centroid and lies outside the physical cross-section.")

if __name__ == '__main__':
    calculate_shear_center()