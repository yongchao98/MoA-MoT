import math

def calculate_shear_center_channel(h, b, t_w, t_f):
    """
    Calculates the shear center location for a symmetric, thin-walled channel section.

    Args:
        h (float): Height of the web (distance between flange centerlines).
        b (float): Width of the flange (from web centerline to flange tip).
        t_w (float): Thickness of the web.
        t_f (float): Thickness of the flange.
    """

    print("--- Shear Center Calculation for a Symmetric Channel Section ---")
    print(f"This demonstrates how the shear center is located outside the section at a distance 'e',")
    print("determined by the section's dimensions, to prevent twisting.")
    print("\nGiven Dimensions:")
    print(f"Web Height (h): {h}")
    print(f"Flange Width (b): {b}")
    print(f"Web Thickness (t_w): {t_w}")
    print(f"Flange Thickness (t_f): {t_f}\n")

    # Step 1: Calculate the moment of inertia (I_x) about the horizontal centroidal axis
    # The formula ignores the flanges' own moment of inertia about their centroids, a common thin-wall approximation.
    I_web = (t_w * h**3) / 12
    # Using Parallel Axis Theorem for the two flanges
    I_flanges = 2 * (b * t_f * (h / 2)**2)
    I_x = I_web + I_flanges

    print("--- Step 1: Calculate Moment of Inertia (I_x) ---")
    print(f"I_x = (t_w * h^3) / 12 + 2 * (b * t_f * (h/2)^2)")
    print(f"I_x = ({t_w} * {h}^3) / 12 + 2 * ({b} * {t_f} * ({h}/2)^2)")
    print(f"I_x = {I_web:.2f} + {I_flanges:.2f} = {I_x:.2f}\n")

    # Step 2: Calculate the shear center offset 'e' from the web centerline
    # This formula balances the external moment (V*e) with the internal moment from shear flow in flanges.
    # Internal Moment M = (V * t_f * h^2 * b^2) / (4 * I_x)
    # V*e = M  => e = M/V
    e = (t_f * h**2 * b**2) / (4 * I_x)

    print("--- Step 2: Calculate Shear Center offset (e) ---")
    print(f"The offset 'e' is found by balancing the external moment with the internal moment from flange shear forces.")
    print(f"e = (t_f * h^2 * b^2) / (4 * I_x)")
    numerator = t_f * h**2 * b**2
    denominator = 4 * I_x
    print(f"e = ({t_f} * {h}^2 * {b}^2) / (4 * {I_x:.2f})")
    print(f"e = {numerator:.2f} / {denominator:.2f}")
    print(f"e = {e:.2f}\n")
    
    print("--- Conclusion ---")
    print(f"The shear center is located {e:.2f} units away from the web's centerline.")
    print("This confirms that the location is outside the section and is a function of its dimensions, which prevents twisting under a shear load.")


# --- Example Usage ---
# Dimensions for an example channel section (e.g., in mm)
web_height = 100.0
flange_width = 50.0
web_thickness = 4.0
flange_thickness = 6.0

calculate_shear_center_channel(web_height, flange_width, web_thickness, flange_thickness)