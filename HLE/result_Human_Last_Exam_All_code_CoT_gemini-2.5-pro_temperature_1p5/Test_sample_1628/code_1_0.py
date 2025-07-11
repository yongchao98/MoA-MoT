import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for an asymmetric channel section.

    The shear center is the point through which an applied shear force
    produces no twisting moment on the section. For a channel section, this
    point is offset from the web.

    This script calculates the location of the horizontal centroidal axis (y_bar),
    the moment of inertia about that axis (I_x), and finally the shear
    center offset (e) from the centerline of the web.
    """
    # --- Define Section Dimensions (in mm) ---
    h = 200.0   # Height between flange centerlines
    t_w = 8.0   # Web thickness
    b_1 = 100.0 # Top flange width
    t_1 = 12.0  # Top flange thickness
    b_2 = 80.0  # Bottom flange width
    t_2 = 10.0  # Bottom flange thickness

    print("--- Section Properties ---")
    print(f"Web Height (h): {h} mm")
    print(f"Web Thickness (t_w): {t_w} mm")
    print(f"Top Flange Width (b_1): {b_1} mm")
    print(f"Top Flange Thickness (t_1): {t_1} mm")
    print(f"Bottom Flange Width (b_2): {b_2} mm")
    print(f"Bottom Flange Thickness (t_2): {t_2} mm")
    print("-" * 28 + "\n")

    # --- Step 1: Calculate Area and Centroid (y_bar) ---
    # The origin (y=0) is at the mid-height of the web.
    A_1 = b_1 * t_1  # Area of top flange
    y_1 = h / 2.0    # Centroid of top flange

    A_2 = b_2 * t_2  # Area of bottom flange
    y_2 = -h / 2.0   # Centroid of bottom flange

    A_w = h * t_w    # Area of web
    y_w = 0.0        # Centroid of web

    A_total = A_1 + A_2 + A_w
    
    # y_bar is the location of the horizontal centroidal axis from the web's mid-height
    y_bar = (A_1 * y_1 + A_2 * y_2 + A_w * y_w) / A_total

    print("--- Centroid Calculation ---")
    print(f"Total Area (A_total): {A_total:.2f} mm^2")
    print(f"Centroid location (y_bar) from web mid-height: {y_bar:.4f} mm")
    print("-" * 28 + "\n")

    # --- Step 2: Calculate Moment of Inertia (I_x) about centroidal axis ---
    # Using the Parallel Axis Theorem: I = I_c + A*d^2
    # I_c for a rectangle is (b*h^3)/12
    I_c1 = (b_1 * t_1**3) / 12.0
    I_1 = I_c1 + A_1 * (y_1 - y_bar)**2

    I_c2 = (b_2 * t_2**3) / 12.0
    I_2 = I_c2 + A_2 * (y_2 - y_bar)**2

    I_cw = (t_w * h**3) / 12.0
    I_w = I_cw + A_w * (y_w - y_bar)**2

    I_x = I_1 + I_2 + I_w

    print("--- Moment of Inertia Calculation (I_x) ---")
    print(f"I_x for Top Flange: {I_1:,.2f} mm^4")
    print(f"I_x for Bottom Flange: {I_2:,.2f} mm^4")
    print(f"I_x for Web: {I_w:,.2f} mm^4")
    print(f"Total Moment of Inertia (I_x): {I_x:,.2f} mm^4")
    print("-" * 28 + "\n")
    
    # --- Step 3: Calculate Shear Center offset (e) from web centerline ---
    # This formula is derived by balancing the external shear moment (V*e)
    # with the internal moment from shear flows in the flanges.
    term1 = (h / 2.0 - y_bar) * t_1 * b_1**2
    term2 = (h / 2.0 + y_bar) * t_2 * b_2**2
    
    e = (h / (4.0 * I_x)) * (term1 + term2)

    print("--- Shear Center Calculation ---")
    print("Equation: e = (h / (4 * I_x)) * [ (h/2 - y_bar)*t_1*b_1^2 + (h/2 + y_bar)*t_2*b_2^2 ]")
    print(f"e = ({h} / (4 * {I_x:.2f})) * [ ({h/2} - {y_bar:.4f})*{t_1}*{b_1**2} + ({h/2} + {y_bar:.4f})*{t_2}*{b_2**2} ]")
    print(f"Shear Center offset (e) from web centerline: {e:.4f} mm")
    print("\nThe shear center is located outside the section, offset from the web by this distance.")

calculate_shear_center_asymmetric_channel()