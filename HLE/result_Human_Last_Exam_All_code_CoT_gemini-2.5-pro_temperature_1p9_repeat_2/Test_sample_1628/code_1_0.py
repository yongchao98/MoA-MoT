import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the approximate location of the shear center for a thin-walled,
    asymmetric channel section.

    The shear center is the point through which applied shear forces do not cause twisting.
    For a channel section, this point lies outside the section.

    This calculation finds the horizontal distance (e_x) of the shear center from the
    centerline of the web. It uses a common thin-wall approximation where the vertical
    centroid is assumed to be at the web's mid-height.
    """
    # --- Define Section Dimensions (in mm) ---
    h = 150.0      # Height of the web
    t_w = 6.0      # Thickness of the web
    b1 = 75.0      # Width of the top flange
    t_f1 = 10.0    # Thickness of the top flange
    b2 = 50.0      # Width of the bottom flange
    t_f2 = 8.0     # Thickness of the bottom flange

    print("--- Asymmetric Channel Section Properties ---")
    print(f"Web Height (h): {h} mm")
    print(f"Web Thickness (t_w): {t_w} mm")
    print(f"Top Flange Width (b1): {b1} mm")
    print(f"Top Flange Thickness (t_f1): {t_f1} mm")
    print(f"Bottom Flange Width (b2): {b2} mm")
    print(f"Bottom Flange Thickness (t_f2): {t_f2} mm")
    print("-" * 43)

    # --- Calculation Steps ---

    # 1. Approximate Moment of Inertia (I_x) about the horizontal axis
    #    at mid-height (h/2). This ignores the flanges' own moment of
    #    inertia about their centroid and assumes t << h.
    #    Formula: I_x = I_web + I_flange1 + I_flange2
    #    where I_flange = A_flange * d^2
    I_web = (t_w * h**3) / 12
    I_flange1 = (b1 * t_f1) * (h / 2)**2
    I_flange2 = (b2 * t_f2) * (h / 2)**2
    I_x_approx = I_web + I_flange1 + I_flange2

    # 2. Calculate the horizontal eccentricity (e_x) of the shear center from the web centerline.
    #    This formula calculates the moment produced by the shear flow in the flanges.
    #    Formula: e_x = (h^2 / (4 * I_x)) * (t_f1 * b1^2 + t_f2 * b2^2)
    numerator = (h**2 / 4) * (t_f1 * b1**2 + t_f2 * b2**2)
    e_x = numerator / I_x_approx
    
    # --- Final Equation and Result ---
    print("Final Calculation for Shear Center Eccentricity (e_x):\n")
    
    # Numerator calculation for the print statement
    num_part1 = h**2 / 4
    num_part2_val_1 = t_f1 * b1**2
    num_part2_val_2 = t_f2 * b2**2
    num_part2 = num_part2_val_1 + num_part2_val_2
    
    print(f"I_x ≈ (t_w*h³/12) + (b1*t_f1*(h/2)²) + (b2*t_f2*(h/2)²)")
    print(f"I_x ≈ ({t_w:.1f}*{h:.1f}³/12) + ({b1:.1f}*{t_f1:.1f}*({h/2:.1f})²) + ({b2:.1f}*{t_f2:.1f}*({h/2:.1f})²)")
    print(f"I_x ≈ {I_web:.2e} + {I_flange1:.2e} + {I_flange2:.2e} = {I_x_approx:.2e} mm⁴\n")

    print(f"e_x = (h² / (4 * I_x)) * (t_f1*b1² + t_f2*b2²)")
    print(f"e_x = ( {h:.1f}² / (4 * {I_x_approx:.2e}) ) * ( {t_f1:.1f}*{b1:.1f}² + {t_f2:.1f}*{b2:.1f}² )")
    print(f"e_x = ( {num_part1:.2e} / {I_x_approx:.2e} ) * ( {num_part2_val_1:.2e} + {num_part2_val_2:.2e} )")
    print(f"e_x = ( {num_part1 / I_x_approx:.4f} ) * ( {num_part2:.2e} )")

    print("\n--- Result ---")
    print(f"The shear center is located at a horizontal distance e_x ≈ {e_x:.2f} mm from the centerline of the web.")
    print("This location is outside the physical cross-section, on the side opposite the flanges.")

if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()