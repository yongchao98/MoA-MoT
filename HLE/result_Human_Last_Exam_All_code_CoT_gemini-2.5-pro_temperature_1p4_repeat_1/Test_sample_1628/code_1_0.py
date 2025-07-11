import math

def calculate_shear_center_asymmetric_channel():
    """
    Calculates the shear center location for an asymmetric channel section.
    The shear center is the point where a shear force can be applied without causing torsion.
    For a channel, this is offset horizontally from the web.
    """
    # --- 1. Define Section Dimensions (in mm) ---
    H = 200.0  # Total height of the section
    b1 = 80.0   # Top flange width
    tf1 = 12.0  # Top flange thickness
    b2 = 100.0  # Bottom flange width
    tf2 = 10.0  # Bottom flange thickness
    tw = 8.0    # Web thickness
    
    print("--- Input Dimensions (mm) ---")
    print(f"Total Height (H): {H}")
    print(f"Top Flange Width (b1): {b1}")
    print(f"Top Flange Thickness (tf1): {tf1}")
    print(f"Bottom Flange Width (b2): {b2}")
    print(f"Bottom Flange Thickness (tf2): {tf2}")
    print(f"Web Thickness (tw): {tw}\n")

    # --- 2. Calculate Centroid (y_c from the bottom outer edge) ---
    # Areas of components
    A1 = b1 * tf1       # Top flange
    A2 = (H - tf1 - tf2) * tw # Web
    A3 = b2 * tf2       # Bottom flange
    A_total = A1 + A2 + A3

    # y-distances to component centroids from the bottom edge
    y1 = H - tf1 / 2.0
    y2 = (H - tf1 - tf2) / 2.0 + tf2
    y3 = tf2 / 2.0
    
    # Vertical centroid y_c
    y_c = (A1 * y1 + A2 * y2 + A3 * y3) / A_total
    
    print("--- Centroid Calculation ---")
    print(f"Area_total = {A1:.2f} + {A2:.2f} + {A3:.2f} = {A_total:.2f} mm^2")
    print(f"y_c = ({A1:.2f}*{y1:.2f} + {A2:.2f}*{y2:.2f} + {A3:.2f}*{y3:.2f}) / {A_total:.2f}")
    print(f"Vertical Centroid (y_c from bottom): {y_c:.2f} mm\n")

    # --- 3. Calculate Moment of Inertia (I_x about the centroidal axis) ---
    # Using Parallel Axis Theorem: I = I_local + A*d^2
    # I_local for rectangles: b*h^3 / 12
    I_local1 = b1 * tf1**3 / 12.0
    I_local2 = tw * (H - tf1 - tf2)**3 / 12.0
    I_local3 = b2 * tf2**3 / 12.0
    
    d1 = y1 - y_c
    d2 = y2 - y_c
    d3 = y3 - y_c
    
    I_x = (I_local1 + A1 * d1**2) + (I_local2 + A2 * d2**2) + (I_local3 + A3 * d3**2)
    
    print("--- Moment of Inertia (I_x) Calculation ---")
    print(f"I_x = ({I_local1:.2f} + {A1:.2f}*{d1:.2f}^2) + ({I_local2:.2f} + {A2:.2f}*{d2:.2f}^2) + ({I_local3:.2f} + {A3:.2f}*{d3:.2f}^2)")
    print(f"Moment of Inertia (I_x): {I_x:.2e} mm^4\n")

    # --- 4. Calculate Shear Center Offset (e from web centerline) ---
    # Formula: e = (1 / (2 * I_x)) * (tf1*h1^2*b1^2 + tf2*h2^2*b2^2)
    # h1 and h2 are distances from the centroidal axis to the flange centerlines
    h1 = H - y_c - tf1 / 2.0 # Note: this is the same as d1
    h2 = y_c - tf2 / 2.0      # Note: this is |d3|
    
    # Calculate e
    numerator = (tf1 * h1**2 * b1**2) + (tf2 * h2**2 * b2**2)
    e = numerator / (2 * I_x)

    print("--- Shear Center Offset (e) Calculation ---")
    print("Formula: e = (tf1*h1^2*b1^2 + tf2*h2^2*b2^2) / (2*I_x)")
    print(f"e = ({tf1:.2f}*{h1:.2f}^2*{b1:.2f}^2 + {tf2:.2f}*{h2:.2f}^2*{b2:.2f}^2) / (2*{I_x:.2e})")
    print(f"e = {numerator:.2e} / {2*I_x:.2e}")
    print(f"Shear Center Offset (e from web centerline): {e:.2f} mm")

if __name__ == '__main__':
    calculate_shear_center_asymmetric_channel()