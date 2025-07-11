import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for an asymmetric channel section.
    The shear center is the point through which shear loads must act to produce
    bending without twisting. For a channel section, it is located outside the
    section, at a distance 'e' from the web.
    """
    # 1. Define section geometry (in mm)
    b1 = 100.0  # Width of top flange
    b2 = 50.0   # Width of bottom flange
    h_total = 200.0 # Total height of the section
    t = 10.0    # Thickness of flanges and web (assumed constant)

    print("--- Geometric Properties ---")
    print(f"Top flange width (b1): {b1} mm")
    print(f"Bottom flange width (b2): {b2} mm")
    print(f"Total height (h_total): {h_total} mm")
    print(f"Thickness (t): {t} mm\n")

    # 2. Calculate centroid location (y_bar) from the bottom edge
    # The section is divided into three rectangles:
    # Part 1: Bottom flange
    area1 = b2 * t
    y1 = t / 2.0
    # Part 2: Web
    h_web = h_total - 2 * t
    area2 = h_web * t
    y2 = t + h_web / 2.0
    # Part 3: Top flange
    area3 = b1 * t
    y3 = h_total - t / 2.0

    total_area = area1 + area2 + area3
    sum_ay = (area1 * y1) + (area2 * y2) + (area3 * y3)
    y_bar = sum_ay / total_area

    print("--- Centroid Calculation ---")
    print(f"Centroid location (y_bar) from bottom edge: {y_bar:.2f} mm\n")

    # 3. Calculate Moment of Inertia (I_x) about the horizontal centroidal axis
    # Using the parallel axis theorem: I_x = sum(I_c + A*d^2)
    # Bottom flange
    I_c1 = (b2 * t**3) / 12.0
    d1 = y_bar - y1
    I_x1 = I_c1 + area1 * d1**2
    # Web
    I_c2 = (t * h_web**3) / 12.0
    d2 = y_bar - y2
    I_x2 = I_c2 + area2 * d2**2
    # Top flange
    I_c3 = (b1 * t**3) / 12.0
    d3 = y3 - y_bar
    I_x3 = I_c3 + area3 * d3**2

    I_x = I_x1 + I_x2 + I_x3

    print("--- Moment of Inertia Calculation ---")
    print(f"Moment of Inertia (I_x) about centroidal axis: {I_x:,.2f} mm^4\n")

    # 4. Calculate shear center offset 'e' from the web centerline
    # Formula: e = (t / (2 * I_x)) * (b1^2 * h1^2 + b2^2 * h2^2)
    # h1 and h2 are distances from the centroidal axis to flange centerlines
    h1 = y3 - y_bar  # same as d3
    h2 = y_bar - y1  # same as d1

    # Print the variables used in the final equation
    print("--- Shear Center Calculation ---")
    print("The final equation for the shear center offset 'e' is:")
    print("e = (t / (2 * I_x)) * (b1^2 * h1^2 + b2^2 * h2^2)\n")
    
    print("Values for the equation:")
    print(f"t = {t}")
    print(f"I_x = {I_x:.2f}")
    print(f"b1 = {b1}")
    print(f"h1 = {h1:.2f} (distance from centroid to top flange centerline)")
    print(f"b2 = {b2}")
    print(f"h2 = {h2:.2f} (distance from centroid to bottom flange centerline)\n")
    
    # Calculate e
    e = (t / (2 * I_x)) * (b1**2 * h1**2 + b2**2 * h2**2)
    
    # Print the final equation with numbers
    print("Final equation with numerical values:")
    final_eq_str = f"e = ({t} / (2 * {I_x:.2f})) * ({b1}**2 * {h1:.2f}**2 + {b2}**2 * {h2:.2f}**2)"
    print(final_eq_str)
    
    print("\n--- Result ---")
    print(f"The shear center is located at a distance e = {e:.2f} mm from the web centerline.")
    print("This location is outside the cross-section.")

if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()