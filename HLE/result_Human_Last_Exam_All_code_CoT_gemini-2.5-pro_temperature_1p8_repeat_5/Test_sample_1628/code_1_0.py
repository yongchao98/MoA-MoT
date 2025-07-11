import math

def calculate_channel_properties(h, tw, b1, tf1, b2, tf2):
    """
    Calculates the centroid and shear center for an asymmetric channel section.

    The origin (0,0) is at the outer corner of the web and the bottom flange's outer face.
    The web lies along the y-axis from y=0 to y=h.

    Args:
        h (float): Height of the web.
        tw (float): Thickness of the web.
        b1 (float): Width of the top flange.
        tf1 (float): Thickness of the top flange.
        b2 (float): Width of the bottom flange.
        tf2 (float): Thickness of the bottom flange.
    """

    print("--- Input Section Dimensions ---")
    print(f"Web Height (h): {h}")
    print(f"Web Thickness (tw): {tw}")
    print(f"Top Flange Width (b1): {b1}")
    print(f"Top Flange Thickness (tf1): {tf1}")
    print(f"Bottom Flange Width (b2): {b2}")
    print(f"Bottom Flange Thickness (tf2): {tf2}\n")

    # Define the three rectangular parts of the section
    # part_1: top flange, part_2: web, part_3: bottom flange
    parts = [
        {'name': 'Top Flange',   'area': b1 * tf1, 'xc': b1 / 2, 'yc': h + tf1 / 2},
        {'name': 'Web',          'area': h * tw,   'xc': tw / 2, 'yc': h / 2},
        {'name': 'Bottom Flange','area': b2 * tf2, 'xc': b2 / 2, 'yc': -tf2 / 2}
    ]

    # --- 1. Calculate Centroid (x_bar, y_bar) ---
    total_area = sum(p['area'] for p in parts)
    sum_ax = sum(p['area'] * p['xc'] for p in parts)
    sum_ay = sum(p['area'] * p['yc'] for p in parts)

    x_bar = sum_ax / total_area
    y_bar = sum_ay / total_area

    print("--- 1. Centroid Calculation ---")
    print(f"Total Area: {total_area:.2f}")
    print(f"Centroid Location (x_bar, y_bar): ({x_bar:.2f}, {y_bar:.2f})\n")

    # --- 2. Calculate Moment of Inertia (Ix) about the centroidal x-axis ---
    ix_total = 0
    # Ix for a rectangle about its own centroid = (b*h^3)/12
    # Parallel Axis Theorem: I = I_c + A*d^2
    ix_total += (b1 * tf1**3 / 12) + parts[0]['area'] * (parts[0]['yc'] - y_bar)**2
    ix_total += (tw * h**3 / 12)   + parts[1]['area'] * (parts[1]['yc'] - y_bar)**2
    ix_total += (b2 * tf2**3 / 12) + parts[2]['area'] * (parts[2]['yc'] - y_bar)**2

    print("--- 2. Moment of Inertia (Ix) ---")
    print(f"Moment of Inertia about Neutral Axis (Ix): {ix_total:.2f}\n")

    # --- 3. Calculate Shear Center Location ---
    # The shear center lies on the neutral axis (y_sc = y_bar) for vertical shear
    y_sc = y_bar

    # Distances from the section's neutral axis to the centroid of each flange
    d1 = abs(parts[0]['yc'] - y_bar) # For top flange
    d2 = abs(parts[2]['yc'] - y_bar) # For bottom flange

    # Eccentricity 'e' is the horizontal offset of the shear center from the web's centerline
    # Formula derived from balancing the external shear moment with the internal moment from flange shear flows
    e = (1 / (2 * ix_total)) * (tf1 * (b1**2) * (d1**2) + tf2 * (b2**2) * (d2**2))

    # The web centerline is at x = tw / 2. Flanges cause twisting "around" the web.
    # The shear center is located on the opposite side of the web to counteract this.
    x_sc = (tw / 2) - e
    
    print("--- 3. Shear Center Calculation ---")
    print("The shear center is the point where a shear force can be applied without causing torsion.")
    print("For a channel section, it is located outside the section, on the opposite side of the web from the flanges.")
    print(f"Distance from top flange centerline to Neutral Axis (d1): {d1:.2f}")
    print(f"Distance from bottom flange centerline to Neutral Axis (d2): {d2:.2f}")
    print(f"Calculated eccentricity from web centerline (e): {e:.2f}")
    print(f"Shear Center Location (x_sc, y_sc): ({x_sc:.2f}, {y_sc:.2f})\n")
    
    print("--- Final Equation for Shear Center Eccentricity ---")
    print("e = (1 / (2 * Ix)) * [tf1 * b1^2 * d1^2 + tf2 * b2^2 * d2^2]")
    print(f"e = (1 / (2 * {ix_total:.2f})) * [{tf1} * {b1}^2 * {d1:.2f}^2 + {tf2} * {b2}^2 * {d2:.2f}^2]")
    print(f"e = {e:.2f}")


if __name__ == '__main__':
    # Define an example asymmetric channel section
    calculate_channel_properties(h=200, tw=10, b1=100, tf1=15, b2=75, tf2=12)