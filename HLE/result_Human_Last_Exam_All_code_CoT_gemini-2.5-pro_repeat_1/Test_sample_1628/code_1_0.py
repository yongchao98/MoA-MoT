import math

def calculate_shear_center_asymmetric_channel(h, b1, b2, t):
    """
    Calculates the shear center for an asymmetric thin-walled channel section.

    The shear center is the point through which shear loads must act to produce
    bending without twisting. For a channel section, this point lies outside
    the section, offset from the web.

    Args:
        h (float): Height of the web (centerline to centerline).
        b1 (float): Width of the top flange.
        b2 (float): Width of the bottom flange.
        t (float): Thickness of the flanges and web (assumed constant).
    """
    print(f"Calculating Shear Center for an Asymmetric Channel with:")
    print(f"  Web height (h) = {h}")
    print(f"  Top flange width (b1) = {b1}")
    print(f"  Bottom flange width (b2) = {b2}")
    print(f"  Thickness (t) = {t}\n")

    # Step 1: Calculate the location of the centroid (y_c) from the centerline of the bottom flange.
    # We use the principle of moments for areas.
    A_top = b1 * t
    A_web = h * t
    A_bot = b2 * t
    A_total = A_top + A_web + A_bot

    # Moment of areas about the bottom flange centerline
    M_areas = (A_top * h) + (A_web * h / 2.0) + (A_bot * 0)
    y_c = M_areas / A_total
    print(f"Step 1: Calculated centroid location (y_c) from bottom flange = {y_c:.4f}\n")

    # Step 2: Calculate the Moment of Inertia (I_x) about the centroidal x-axis.
    # Using the parallel axis theorem: I = I_c + A*d^2
    # The term (1/12)bt^3 for flanges is negligible for thin sections and ignored.
    d_top = h - y_c
    d_web = h / 2.0 - y_c
    d_bot = -y_c

    I_top = A_top * (d_top**2)
    I_web = (1/12.0) * t * (h**3) + A_web * (d_web**2)
    I_bot = A_bot * (d_bot**2)
    I_x = I_top + I_web + I_bot
    print(f"Step 2: Calculated Moment of Inertia (I_x) = {I_x:.4f}\n")

    # Step 3: Calculate the shear center offset (e) from the web centerline.
    # This is derived by balancing the external moment (V*e) with the internal
    # moment from the shear flow in the flanges.
    # Internal Moment M_int = (V*t / (2*I_x)) * [ b1^2 * d_top^2 + b2^2 * d_bot^2 ]
    # e = M_int / V
    e = (t / (2.0 * I_x)) * ( (b1**2) * (d_top**2) + (b2**2) * (d_bot**2) )
    print(f"Step 3: Calculated Shear Center offset (e) from web = {e:.4f}\n")
    
    print("--- Conclusion ---")
    print(f"The shear center is located at a distance 'e' = {e:.4f} outside the web.")
    print("This confirms that its location is not at the centroid, but is a specific distance")
    print("calculated based on the channel's dimensions (flange widths, web height, thickness).")
    print("Applying a shear force at this point ensures the channel bends without twisting.")


# --- Example Calculation ---
# Define the dimensions of a specific asymmetric channel section
web_height = 200.0  # mm
top_flange_width = 80.0   # mm
bottom_flange_width = 120.0 # mm
thickness = 10.0      # mm

# Run the calculation and print the results
calculate_shear_center_asymmetric_channel(web_height, top_flange_width, bottom_flange_width, thickness)