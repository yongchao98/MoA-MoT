import math

def calculate_shear_center_symmetric_channel(h, b, t_w, t_f):
    """
    Calculates the shear center location for a symmetric thin-walled channel section.

    The shear center is the point through which a shear force can be applied
    without causing torsion (twisting). For a symmetric channel (C-section),
    the shear center lies on the axis of symmetry, outside the section, a
    distance 'e' from the web centerline.

    For an asymmetric channel, the principle is the same but the calculation is
    more complex as the centroid must be found first, and the contributions
    from both unequal flanges must be accounted for.

    Args:
        h (float): Height of the web (distance between flange centerlines).
        b (float): Width of the flange (from web centerline to tip).
        t_w (float): Thickness of the web.
        t_f (float): Thickness of the flanges.
    """
    print("--- Input Section Dimensions ---")
    print(f"Web height (h): {h}")
    print(f"Flange width (b): {b}")
    print(f"Web thickness (t_w): {t_w}")
    print(f"Flange thickness (t_f): {t_f}\n")

    # --- Step 1: Calculate the moment of inertia (I_x) about the x-axis ---
    # Using the parallel axis theorem, approximating for a thin-walled section.
    # I_x = I_web + I_flanges
    # I_x for web (passes through centroid) = (t_w * h^3) / 12
    # I_x for one flange = I_c + A*d^2 ~= (b*t_f^3)/12 + (b*t_f) * (h/2)^2
    # We use the common engineering approximation which neglects the I_c of the flange.
    i_web = (t_w * (h**3)) / 12
    # For both flanges
    i_flanges = 2 * (b * t_f * ((h/2)**2))
    i_x = i_web + i_flanges

    print("--- Calculation Steps ---")
    print(f"1. Moment of Inertia (I_x) = (t_w * h^3)/12 + 2 * b * t_f * (h/2)^2")
    print(f"   I_x = ({t_w} * {h}^3)/12 + 2 * {b} * {t_f} * ({h}/2)^2")
    print(f"   I_x = {i_web:.2f} + {i_flanges:.2f} = {i_x:.2f}\n")


    # --- Step 2: Calculate the shear center offset (e) from the web centerline ---
    # This formula is derived from balancing the external shear moment with the
    # internal moment from the shear flow in the flanges.
    # e = (b^2 * h^2 * t_f) / (4 * I_x)
    numerator = (b**2) * (h**2) * t_f
    denominator = 4 * i_x
    e = numerator / denominator
    
    print(f"2. Shear Center Offset (e) = (b^2 * h^2 * t_f) / (4 * I_x)")
    print(f"   e = ({b}^2 * {h}^2 * {t_f}) / (4 * {i_x:.2f})")
    print(f"   e = {numerator:.2f} / {denominator:.2f}\n")


    print("--- Result ---")
    print(f"The shear center is located at a distance e = {e:.4f} from the centerline of the web.")
    print("This confirms the location is outside the section at a calculated distance based on dimensions.")

# --- Example Usage ---
# Define the dimensions for a symmetric channel section (e.g., in mm)
web_height = 150.0
flange_width = 75.0
web_thickness = 5.0
flange_thickness = 8.0

calculate_shear_center_symmetric_channel(web_height, flange_width, web_thickness, flange_thickness)