import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for an asymmetric channel section.

    The shear center is the point where a shear force can be applied
    without causing torsion. For a channel section, it is located
    outside the section, offset from the web.

    This calculation assumes:
    - The section is thin-walled.
    - The vertical shear force Vy is applied.
    - The horizontal centroidal axis is at the mid-height of the web (a common simplification).
    """

    # --- 1. Define Geometric Properties (in mm) ---
    h = 150.0   # Height of the web
    t_w = 6.0   # Thickness of the web
    b1 = 80.0   # Width of the top flange
    b2 = 50.0   # Width of the bottom flange
    t_f = 10.0  # Thickness of the flanges (assumed uniform)

    print("Geometric Properties of the Asymmetric Channel Section:")
    print(f"Web height (h): {h} mm")
    print(f"Web thickness (t_w): {t_w} mm")
    print(f"Top flange width (b1): {b1} mm")
    print(f"Bottom flange width (b2): {b2} mm")
    print(f"Flange thickness (t_f): {t_f} mm\n")

    # --- 2. Calculate Moment of Inertia (I_x) about the horizontal centroidal axis ---
    # Simplified formula neglecting the flanges' own inertia about their centroid (I = I_web + A_flanges * d^2)
    I_web = (t_w * h**3) / 12
    I_flanges = (b1 * t_f * (h/2)**2) + (b2 * t_f * (h/2)**2)
    I_x = I_web + I_flanges

    print("--- Calculation Steps ---")
    print("1. Calculate Moment of Inertia (I_x):")
    print(f"I_x = (t_w * h^3)/12 + (b1*t_f + b2*t_f) * (h/2)^2")
    print(f"I_x = ({t_w} * {h}^3)/12 + ({b1}*{t_f} + {b2}*{t_f}) * ({h}/2)^2")
    print(f"I_x = {I_web:.2f} + {I_flanges:.2f} = {I_x:.2f} mm^4\n")

    # --- 3. Calculate Shear Center offset (e) from the web ---
    # Formula derived from balancing moments: V*e = Moment from flange shear flows
    numerator = t_f * h**2 * (b1**2 + b2**2)
    denominator = 4 * I_x
    e = numerator / denominator
    
    print("2. Calculate Shear Center offset (e) from the web:")
    print(f"e = (t_f * h^2 * (b1^2 + b2^2)) / (4 * I_x)")
    # Print the equation with numbers
    print(f"e = ({t_f} * {h}^2 * ({b1}^2 + {b2}^2)) / (4 * {I_x:.2f})")
    print(f"e = ({t_f} * {h**2} * ({b1**2} + {b2**2})) / ({4 * I_x:.2f})")
    print(f"e = ({t_f * h**2 * (b1**2 + b2**2):.2f}) / ({4 * I_x:.2f})")
    print(f"\nFinal Result:")
    print(f"The shear center is located at a horizontal distance e = {e:.2f} mm from the web.")


if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()