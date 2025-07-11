import math

def calculate_shear_center_asymmetric_channel(h, tw, b1, t1, b2, t2):
    """
    Calculates the shear center location for an asymmetric channel section.

    The shear center is the point where a shear load can be applied without
    causing torsion. For a channel, it is located horizontally outside the
    section, offset from the web.

    Args:
        h (float): Height of the web (distance between flange centerlines).
        tw (float): Thickness of the web.
        b1 (float): Width of the top flange.
        t1 (float): Thickness of the top flange.
        b2 (float): Width of the bottom flange.
        t2 (float): Thickness of the bottom flange.
    """

    print("1. Define Section Properties:")
    print(f"  Web Height (h): {h}")
    print(f"  Web Thickness (tw): {tw}")
    print(f"  Top Flange (b1 x t1): {b1} x {t1}")
    print(f"  Bottom Flange (b2 x t2): {b2} x {t2}\n")

    # --- Step 1: Calculate vertical centroid (y_c) ---
    # Origin is at the mid-height of the web. Positive y is up.
    A_web = h * tw
    y_web_c = 0.0

    A_top_flange = b1 * t1
    y_top_flange_c = h / 2.0

    A_bottom_flange = b2 * t2
    y_bottom_flange_c = -h / 2.0

    A_total = A_web + A_top_flange + A_bottom_flange
    
    # Sum of moments of area about the web's mid-height
    sum_Ay = (A_web * y_web_c) + (A_top_flange * y_top_flange_c) + (A_bottom_flange * y_bottom_flange_c)
    
    y_c = sum_Ay / A_total

    print(f"2. Calculate Vertical Centroid (y_c) from web mid-height:")
    print(f"  y_c = (A_top*y_top + A_web*y_web + A_bot*y_bot) / A_total")
    print(f"  y_c = ({A_top_flange:.2f}*{y_top_flange_c:.2f} + {A_web:.2f}*{y_web_c:.2f} + {A_bottom_flange:.2f}*{y_bottom_flange_c:.2f}) / {A_total:.2f}")
    print(f"  y_c = {y_c:.4f}\n")


    # --- Step 2: Calculate Moment of Inertia (I_x) about centroidal axis ---
    # Using the parallel axis theorem: I = I_c + A*d^2
    # For thin rectangles, I_c about axis parallel to width is negligible (e.g., b*t^3/12 is small)
    
    # Distance from overall centroid to each component's centroid
    d_web = y_web_c - y_c
    d_top = y_top_flange_c - y_c
    d_bot = y_bottom_flange_c - y_c

    # I_x for web
    I_x_web = (tw * h**3) / 12 + A_web * d_web**2
    # I_x for top flange (ignoring b1*t1^3/12 term)
    I_x_top = A_top_flange * d_top**2
    # I_x for bottom flange (ignoring b2*t2^3/12 term)
    I_x_bot = A_bottom_flange * d_bot**2

    I_x = I_x_web + I_x_top + I_x_bot

    print("3. Calculate Moment of Inertia (I_x) about Centroidal Axis:")
    print(f"  I_x = I_web + I_top_flange + I_bottom_flange")
    print(f"  I_x = ({I_x_web:.2f}) + ({I_x_top:.2f}) + ({I_x_bot:.2f})")
    print(f"  I_x = {I_x:.4f}\n")


    # --- Step 3: Calculate Shear Center Eccentricity (e) ---
    # e is the horizontal distance from the web centerline to the shear center.
    # The formula is derived by balancing the external shear moment (V*e)
    # with the internal moment from shear flows in the flanges.
    # e = (Moment_from_flanges) / V
    # Moment_from_flanges = F_top * d_top + F_bot * d_bot
    
    # The formula simplifies to:
    # e = (h / (2 * I_x)) * [ (t1*b1^2/2)*d_top + (t2*b2^2/2)*d_bot ] -> My previous derivation was slightly off
    # A more robust formula for horizontal eccentricity from the web centerline is:
    # e = (1 / I_x) * [ (t1 * b1**2 * h * d_top)/4 + (t2 * b2**2 * h * abs(d_bot))/4 ] -> No
    # The moment from top flange flow is F_top*(h/2) and bottom is F_bot*(h/2)
    # The force F_top is proportional to (b1^2/2)*t1*d_top. Force F_bot to (b2^2/2)*t2*abs(d_bot).
    
    term_top = (t1 * b1**2 / 2) * d_top
    term_bot = (t2 * b2**2 / 2) * abs(d_bot)
    # We take absolute value of d_bot because the moment contribution adds up.
    
    # This formula appears in some derivations, but let's use the most direct one: sum of moments from flange shear forces.
    # Moment about centroid M_c = F_top*d_top + F_bot_mag*d_bot
    moment_from_flanges_normalized = (t1 * d_top**2 * b1**2 / 2) + (t2 * d_bot**2 * b2**2 / 2)
    e = (1 / I_x) * moment_from_flanges_normalized

    print("4. Calculate Shear Center Horizontal Eccentricity (e) from Web Centerline:")
    print("  e = (1 / I_x) * [ (t1*d_top^2*b1^2 / 2) + (t2*d_bot^2*b2^2 / 2) ]")
    print(f"  e = (1 / {I_x:.2f}) * [ ({t1}*{d_top:.2f}^2*{b1}^2 / 2) + ({t2}*{d_bot:.2f}^2*{b2}^2 / 2) ]")
    
    val1 = (t1 * d_top**2 * b1**2 / 2)
    val2 = (t2 * d_bot**2 * b2**2 / 2)
    print(f"  e = (1 / {I_x:.2f}) * [ {val1:.2f} + {val2:.2f} ]")
    print(f"  e = {e:.4f}\n")
    
    print("--- Final Location ---")
    print(f"The shear center is located {e:.4f} units horizontally from the web's centerline, outside the section.")


if __name__ == '__main__':
    # Example asymmetric channel section properties
    calculate_shear_center_asymmetric_channel(
        h=100.0,   # web height
        tw=6.0,    # web thickness
        b1=80.0,   # top flange width
        t1=8.0,    # top flange thickness
        b2=40.0,   # bottom flange width
        t2=8.0     # bottom flange thickness
    )