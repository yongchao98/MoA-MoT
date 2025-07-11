import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for a thin-walled asymmetric channel section.
    """
    # --- 1. Define Section Geometry (units in mm) ---
    # h is the clear height of the web, or distance between flange centerlines for simplicity
    h = 150.0  # Height of the web
    b1 = 75.0   # Width of the top flange
    b2 = 50.0   # Width of the bottom flange
    t = 6.0     # Thickness (assumed uniform)

    print("--- Geometric Properties (mm) ---")
    print(f"Web height (h): {h}")
    print(f"Top flange width (b1): {b1}")
    print(f"Bottom flange width (b2): {b2}")
    print(f"Thickness (t): {t}\n")

    # --- 2. Calculate Centroid Location (y_bar) ---
    # y_bar is the vertical distance of the centroid from the centerline of the bottom flange.
    # Area of each component
    area_top_flange = b1 * t
    area_bottom_flange = b2 * t
    area_web = h * t

    # Total Area
    total_area = area_top_flange + area_bottom_flange + area_web

    # Moment of area about the bottom flange's centerline
    moment_area = (area_top_flange * h) + (area_web * h / 2.0)

    # Centroid location
    y_bar = moment_area / total_area

    print("--- Centroid Calculation ---")
    print(f"y_bar (from bottom flange) = ((b1*t*h) + (h*t*h/2)) / (b1*t + b2*t + h*t)")
    print(f"y_bar = (({b1}*{t}*{h}) + ({h}*{t}*{h}/2)) / (({b1}*{t}) + ({b2}*{t}) + ({h}*{t}))")
    print(f"y_bar = {y_bar:.2f} mm\n")
    
    # Vertical distances from the centroid to the flange centerlines
    y_top = h - y_bar
    y_bot = -y_bar # It is below the centroid, so its coordinate is negative
    
    print("--- Distances from Centroid to Flanges ---")
    print(f"Distance to top flange (y_top): {y_top:.2f} mm")
    print(f"Distance to bottom flange (y_bot): {y_bot:.2f} mm\n")


    # --- 3. Calculate Second Moment of Area (I_x) about the centroidal axis ---
    # Using the parallel axis theorem: I_x = Σ(I_c + A*d^2)
    # I_c for flanges (t*b^3/12) is negligible for thin sections rotating about x-axis.
    # I_c for web = t*h^3/12
    
    # Top flange
    I_x_top = (area_top_flange * y_top**2)
    # Bottom flange
    I_x_bot = (area_bottom_flange * y_bot**2)
    # Web
    I_c_web = (t * h**3) / 12.0
    d_web = (h / 2.0) - y_bar
    I_x_web = I_c_web + (area_web * d_web**2)

    I_x = I_x_top + I_x_bot + I_x_web
    
    print("--- Second Moment of Area Calculation (I_x) ---")
    print(f"I_x = I_top_flange + I_bottom_flange + I_web")
    print(f"I_x = (A_top * y_top^2) + (A_bot * y_bot^2) + (t*h^3/12 + A_web * d_web^2)")
    print(f"I_x = (({b1}*{t})*{y_top:.2f}^2) + (({b2}*{t})*{y_bot:.2f}^2) + ({t}*{h}^3/12 + ({h}*{t})*{d_web:.2f}^2)")
    print(f"I_x = {I_x_top:.0f} + {I_x_bot:.0f} + {I_x_web:.0f}")
    print(f"I_x = {I_x:.2f} mm^4\n")


    # --- 4. Calculate Shear Center offset (e) from web centerline ---
    # Formula: e = (t / (2 * I_x)) * (b1^2 * y_top^2 + b2^2 * y_bot^2)
    term1 = b1**2 * y_top**2
    term2 = b2**2 * y_bot**2 # y_bot is squared so the sign doesn't matter
    
    e = (t / (2 * I_x)) * (term1 + term2)

    print("--- Shear Center Offset Calculation (e) ---")
    print(f"The shear center is offset horizontally from the web by a distance 'e'.")
    print(f"e = (t / (2 * I_x)) * (b1^2 * y_top^2 + b2^2 * y_bot^2)")
    print(f"e = ({t} / (2 * {I_x:.2f})) * ({b1}^2 * {y_top:.2f}^2 + {b2}^2 * {y_bot:.2f}^2)")
    print(f"e = ({t / (2 * I_x):.8f}) * ({term1:.2f} + {term2:.2f})")
    print(f"Final Shear Center offset e = {e:.2f} mm")

if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()