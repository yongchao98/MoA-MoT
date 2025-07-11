import math

def calculate_shear_center_symmetric_channel():
    """
    Calculates the shear center location for a SYMMETRIC channel section to demonstrate the principle.
    The shear center is located outside the section at a calculated distance,
    ensuring that an applied shear force does not cause twisting.
    For an asymmetric channel, the principle is the same, but the formulas are more complex.
    """
    # Define the dimensions of the symmetric channel section (in mm)
    h = 150.0  # Overall height of the section
    b = 75.0   # Width of the flanges
    t_w = 6.0    # Thickness of the web
    t_f = 10.0   # Thickness of the flanges

    print("This script calculates the shear center for a symmetric channel section to demonstrate the underlying principle.")
    print("-" * 60)
    print(f"Given Section Dimensions:")
    print(f"  Web Height (h): {h} mm")
    print(f"  Flange Width (b): {b} mm")
    print(f"  Web Thickness (t_w): {t_w} mm")
    print(f"  Flange Thickness (t_f): {t_f} mm")
    print("-" * 60)

    # Use thin-walled approximation where h_centerline is the distance between flange centerlines.
    h_centerline = h - t_f

    # --- Step 1: Calculate Moment of Inertia (I_x) about the strong axis ---
    # I_x = I_web + I_flanges
    # I_flanges uses the parallel axis theorem: I = I_centroidal + A*d^2
    # The centroidal moment of inertia of the thin flange (b*t_f^3/12) is negligible.
    I_x_web = (t_w * h_centerline**3) / 12
    I_x_flanges = 2 * (b * t_f) * (h_centerline / 2)**2
    I_x = I_x_web + I_x_flanges
    
    print("Step 1: Calculate Moment of Inertia (I_x)")
    print("Formula (thin-wall approx.): I_x = (t_w * (h-t_f)^3)/12 + 2 * (b * t_f) * ((h-t_f)/2)^2")
    print(f"Calculation: I_x = ({t_w} * ({h}-{t_f})^3)/12 + 2 * ({b} * {t_f}) * (({h}-{t_f})/2)^2")
    print(f"Result: I_x = {I_x:.2f} mm^4")
    print("-" * 60)

    # --- Step 2: Calculate Shear Center offset 'e' ---
    # The formula is derived by balancing the external moment (V*e) with the internal
    # moment from the shear flow in the flanges.
    # e = (t_f * h_centerline^2 * b^2) / (4 * I_x)
    e = (t_f * h_centerline**2 * b**2) / (4 * I_x)

    print("Step 2: Calculate Shear Center offset (e) from the web centerline")
    print("Formula: e = (t_f * (h-t_f)^2 * b^2) / (4 * I_x)")
    # Final equation with numbers plugged in
    print(f"Calculation: e = ({t_f} * ({h}-{t_f})^2 * {b}^2) / (4 * {I_x:.2f})")
    print(f"Result: e = {e:.2f} mm")
    print("-" * 60)

    print(f"Final Answer:")
    print(f"The shear center is located outside the section at a distance 'e' = {e:.2f} mm from the web.")
    print("This confirms the location is a calculated distance based on the section's dimensions.")

if __name__ == '__main__':
    calculate_shear_center_symmetric_channel()