import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center for a given asymmetric channel section.
    """
    # --- Step 0: Define Section Dimensions (in mm) ---
    # b1, t1: Top flange width and thickness
    # b2, t2: Bottom flange width and thickness
    # h: Total height of the section
    # tw: Web thickness
    
    b1 = 100.0
    t1 = 15.0
    b2 = 75.0
    t2 = 10.0
    h = 250.0
    tw = 8.0

    print("--- Calculating Shear Center for an Asymmetric Channel Section ---")
    print(f"Dimensions: h={h}, b1={b1}, t1={t1}, b2={b2}, t2={t2}, tw={tw}\n")

    # --- Step 1: Calculate Vertical Centroid (y_c) ---
    # We will set the origin (y=0) at the bottommost edge of the section.
    h_web_clear = h - t1 - t2

    # Properties of the top flange
    area_top = b1 * t1
    yc_top = h - t1 / 2.0

    # Properties of the web
    area_web = h_web_clear * tw
    yc_web = t2 + h_web_clear / 2.0

    # Properties of the bottom flange
    area_bot = b2 * t2
    yc_bot = t2 / 2.0
    
    total_area = area_top + area_web + area_bot
    y_centroid = (area_top * yc_top + area_web * yc_web + area_bot * yc_bot) / total_area

    print(f"1. Vertical Centroid Location (from bottom edge): {y_centroid:.2f} mm")

    # --- Step 2: Calculate Moment of Inertia (Ix) about the Centroidal Axis ---
    # Using the Parallel Axis Theorem: I = I_c + A*d^2
    
    # Moment of inertia for the top flange
    I_c_top = (b1 * t1**3) / 12.0
    d_top = yc_top - y_centroid
    Ix_top = I_c_top + area_top * d_top**2

    # Moment of inertia for the web
    I_c_web = (tw * h_web_clear**3) / 12.0
    d_web = yc_web - y_centroid
    Ix_web = I_c_web + area_web * d_web**2
    
    # Moment of inertia for the bottom flange
    I_c_bot = (b2 * t2**3) / 12.0
    d_bot = yc_bot - y_centroid
    Ix_bot = I_c_bot + area_bot * d_bot**2

    Ix_total = Ix_top + Ix_web + Ix_bot
    print(f"2. Moment of Inertia (Ix) about Centroidal Axis: {Ix_total:,.2f} mm^4")

    # --- Step 3: Calculate Shear Center Offset (e) from Web Centerline ---
    # h1 and h2 are distances from the section centroid to the flange centroids
    h1 = abs(d_top)
    h2 = abs(d_bot)
    
    numerator = (t1 * b1**2 * h1**2) + (t2 * b2**2 * h2**2)
    denominator = 2 * Ix_total
    e = numerator / denominator

    print("\n3. Shear Center Offset Calculation (e)")
    print("   Formula: e = (t1*b1^2*h1^2 + t2*b2^2*h2^2) / (2 * Ix)")
    print(f"   Substituting values:")
    print(f"   e = ({t1}*{b1}^2*{h1:.2f}^2 + {t2}*{b2}^2*{h2:.2f}^2) / (2 * {Ix_total:,.2f})")
    print(f"   e = {numerator:,.2f} / {denominator:,.2f}")
    print(f"   e = {e:.2f} mm")
    
    print("\n---")
    print("Conclusion:")
    print("For an asymmetric channel section, the shear center is located outside the cross-section,")
    print("offset from the web by a calculated distance to ensure that applied shear forces do not cause twisting.")
    print(f"For the given dimensions, this offset 'e' is {e:.2f} mm from the web's centerline.")


if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()