def calculate_shear_center_asymmetric_channel():
    """
    Calculates and explains the shear center location for an asymmetric channel section.
    """
    # --- Section Dimensions (user can modify these values) ---
    # All dimensions are in millimeters (mm)
    h = 200.0   # Height of the web
    b1 = 100.0  # Width of the top flange
    b2 = 50.0   # Width of the bottom flange
    t = 10.0    # Uniform thickness of web and flanges

    print("--- Asymmetric Channel Section Properties ---")
    print(f"Web Height (h): {h} mm")
    print(f"Top Flange Width (b1): {b1} mm")
    print(f"Bottom Flange Width (b2): {b2} mm")
    print(f"Uniform Thickness (t): {t} mm\n")

    # --- Step 1: Calculate Centroid (y_bar) ---
    # y_bar is the distance from the centerline of the bottom flange to the centroid.
    area1 = b1 * t  # Top flange area
    y1 = h          # y-distance to top flange centerline
    area2 = h * t   # Web area
    y2 = h / 2      # y-distance to web centerline
    area3 = b2 * t  # Bottom flange area
    y3 = 0.0        # y-distance to bottom flange centerline

    total_area = area1 + area2 + area3
    y_bar = (area1 * y1 + area2 * y2 + area3 * y3) / total_area

    print("--- Calculation Steps ---")
    print(f"1. Vertical Centroid (y_bar) from bottom flange centerline: {y_bar:.2f} mm")

    # --- Step 2: Calculate Moment of Inertia (I_x) ---
    # Using the Parallel Axis Theorem: I_x = sum(I_c + A*d^2) for each part.
    # I_c for thin flanges (b*t^3/12) is negligible.
    
    # For Top Flange (1)
    d1 = h - y_bar
    I_x1 = area1 * d1**2
    
    # For Web (2)
    I_c2 = (t * h**3) / 12
    d2 = (h / 2) - y_bar
    I_x2 = I_c2 + area2 * d2**2
    
    # For Bottom Flange (3)
    d3 = 0 - y_bar
    I_x3 = area3 * d3**2

    I_x = I_x1 + I_x2 + I_x3
    print(f"2. Moment of Inertia about horizontal centroidal axis (I_x): {I_x:.2f} mm^4")

    # --- Step 3: Calculate Shear Center offset (e) from the web centerline ---
    h1_dist = d1     # Distance from centroid to top flange
    h2_dist = y_bar  # Distance from centroid to bottom flange

    # The formula balances the external moment V*e with the internal moment from flange shear flows.
    e = (h * t) / (4 * I_x) * (h1_dist * b1**2 + h2_dist * b2**2)
    
    print("\n3. Shear Center Offset (e) Calculation:")
    print("   Formula: e = (h * t) / (4 * I_x) * (h1 * b1^2 + h2 * b2^2)")
    print("   Substituting values:")
    print(f"   e = ({h:.1f} * {t:.1f}) / (4 * {I_x:.2f}) * ({h1_dist:.2f} * {b1:.1f}^2 + {h2_dist:.2f} * {b2:.1f}^2)")
    
    term1 = (h * t) / (4 * I_x)
    term2_part1 = h1_dist * b1**2
    term2_part2 = h2_dist * b2**2
    term2 = term2_part1 + term2_part2
    
    print(f"   e = {term1:.8f} * ({term2_part1:.2f} + {term2_part2:.2f})")
    print(f"   e = {term1:.8f} * {term2:.2f}")

    print("\n--- Final Result ---")
    print(f"The shear center is located at a distance 'e' = {e:.2f} mm")
    print("This distance is measured horizontally from the centerline of the web, on the side opposite the flanges.")

# Run the calculation
calculate_shear_center_asymmetric_channel()