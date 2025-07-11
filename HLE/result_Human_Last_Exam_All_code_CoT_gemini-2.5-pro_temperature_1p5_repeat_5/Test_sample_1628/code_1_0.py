import math

def calculate_shear_center_channel():
    """
    Calculates the shear center for a singly-symmetric channel section
    to illustrate its location relative to the cross-section.
    
    The principle demonstrated holds true for asymmetric channels as well:
    the shear center is located outside the section.
    """
    
    # --- Define Section Properties (in mm) ---
    # We will use a standard channel section for this example.
    # h is the height of the web between the centerlines of the two flanges.
    # b is the width of the flange, measured from the centerline of the web.
    h = 200.0   # Web height between flange centerlines
    b = 50.0    # Flange width from web centerline
    t_w = 8.0   # Web thickness
    t_f = 12.0  # Flange thickness
    
    print("--- Channel Section Properties ---")
    print(f"Web Height (h): {h} mm")
    print(f"Flange Width (b): {b} mm")
    print(f"Web Thickness (t_w): {t_w} mm")
    print(f"Flange Thickness (t_f): {t_f} mm")
    print("-" * 34)
    
    # --- Step 1: Calculate Moment of Inertia (Ix) about the horizontal centroidal axis ---
    # The formula ignores the flanges' own moment of inertia about their centroid,
    # which is a standard approximation for thin-walled sections.
    # Ix = (Ix of web) + (Ix of flanges using parallel axis theorem)
    
    i_x_web = (t_w * h**3) / 12
    i_x_flanges = 2 * (b * t_f) * (h / 2)**2
    i_x_total = i_x_web + i_x_flanges
    
    print(f"Moment of Inertia (Ix) = {i_x_total:,.2f} mm^4")
    
    # --- Step 2: Calculate Shear Center location 'e' from the web centerline ---
    # This formula calculates the offset required to balance the twisting moment
    # produced by the shear flow in the flanges.
    
    e = (t_f * h**2 * b**2) / (4 * i_x_total)
    
    print("\n--- Shear Center Calculation ---")
    print("Formula: e = (t_f * h^2 * b^2) / (4 * Ix)")
    
    # Print the equation with the numbers plugged in
    print(f"Equation: e = ({t_f} * {h}^2 * {b}^2) / (4 * {i_x_total:,.2f})")
    
    print(f"\nResult: Shear Center offset (e) = {e:.2f} mm")
    
    print("\n--- Conclusion ---")
    print(f"The shear center is located {e:.2f} mm away from the web's centerline.")
    print("This location is outside the physical cross-section.")
    print("When a shear force is applied at this point, it will not cause the channel to twist.")
    print("\nThis confirms that for a channel section, the shear center is located outside the cross-section, offset from the centroid by a distance determined by the channel’s flange widths and web height, ensuring that applied shear forces do not cause twisting.")

# Run the calculation and print the results
calculate_shear_center_channel()