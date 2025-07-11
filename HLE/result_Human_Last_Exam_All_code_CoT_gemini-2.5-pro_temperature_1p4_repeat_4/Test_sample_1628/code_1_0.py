import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for an asymmetric channel section
    using thin-wall approximations.
    """
    # --- Define Section Dimensions (in mm) ---
    b1 = 100.0  # Top flange width
    b2 = 50.0   # Bottom flange width
    h = 200.0   # Height between flange centerlines (web height)
    tf = 10.0   # Flange thickness (assumed constant)
    tw = 8.0    # Web thickness

    print("This script calculates the shear center for an asymmetric channel section.")
    print("The shear center is the point where a shear force can be applied without causing twisting.")
    print("For a channel section, it is located outside the section, away from the web.")
    print("\n--- Given Section Dimensions ---")
    print(f"Top Flange Width (b1): {b1} mm")
    print(f"Bottom Flange Width (b2): {b2} mm")
    print(f"Web Height (h): {h} mm")
    print(f"Flange Thickness (tf): {tf} mm")
    print(f"Web Thickness (tw): {tw} mm")

    # --- Calculation Steps ---
    # 1. Calculate the approximate Moment of Inertia (Ix) about the horizontal
    #    axis passing through the center of the web. This is a common simplification
    #    for thin-walled sections.
    #    Ix = I_web + I_top_flange + I_bottom_flange
    #    I_web ~ (tw * h^3) / 12
    #    I_flange ~ I_flange_centroidal + Area * distance^2
    #    I_flange_centroidal is negligible for thin flanges.
    #    Area_top = b1 * tf, distance_top = h / 2
    #    Area_bottom = b2 * tf, distance_bottom = h / 2
    Ix_approx = (tw * h**3) / 12 + (b1 * tf * (h/2)**2) + (b2 * tf * (h/2)**2)
    Ix_approx = (tw * h**3) / 12 + (tf * h**2 / 4) * (b1 + b2)

    # 2. Calculate the distance 'e' of the shear center from the web's centerline.
    #    This is derived by balancing the moment from the external shear force (V * e)
    #    with the internal moment from the shear flow in the flanges.
    #    Internal Moment = (Force_top_flange + Force_bottom_flange) * (h / 2)
    #    A simplified formula for 'e' is:
    #    e = (tf * h^2 * (b1^2 + b2^2)) / (8 * Ix)
    numerator_e = tf * h**2 * (b1**2 + b2**2)
    e = numerator_e / (8 * Ix_approx)

    print("\n--- Calculation Results ---")
    print(f"1. Approximate Moment of Inertia (Ix) about web's horizontal centerline: {Ix_approx:,.2f} mm^4")
    print(f"2. The horizontal distance 'e' of the shear center from the web centerline is calculated.")
    print(f"   e = (tf * h^2 * (b1^2 + b2^2)) / (8 * Ix)")
    print(f"   e = ({tf:.1f} * {h:.1f}^2 * ({b1:.1f}^2 + {b2:.1f}^2)) / (8 * {Ix_approx:,.2f})")
    print(f"\nFinal Answer: The shear center is located {e:.2f} mm from the web's centerline, outside the section.")
    
    print("\n-------------------------------")
    print("In summary:")
    print("For an asymmetric channel section under pure torsion, the shear center is located outside the cross-section, offset from the centroid by a distance determined by the channel’s flange widths and web height, ensuring that applied shear forces do not cause twisting.")


calculate_shear_center_asymmetric_channel()