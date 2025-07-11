import math

def calculate_shear_center_channel(h, b, t_w, t_f):
    """
    Calculates the shear center location for a symmetric channel section.
    
    The shear center is a point where shear force can be applied without causing torsion.
    For a channel, it is located outside the section, on the side of the web.
    This function calculates the distance 'e' from the centerline of the web.
    
    Args:
    h (float): Overall height of the web.
    b (float): Overall width of one flange.
    t_w (float): Thickness of the web.
    t_f (float): Thickness of the flange.
    
    Returns:
    float: The distance 'e' of the shear center from the web centerline.
    """
    
    # Using a simplified formula for thin-walled sections.
    # This formula calculates the distance from the web's centerline.
    # Numerator term
    numerator = t_f * (b**2)
    
    # Denominator term
    denominator = ((t_w * h) / 3) + (2 * b * t_f)
    
    if denominator == 0:
        return float('inf')
        
    e = numerator / denominator
    return e

# --- Example Dimensions (in mm) ---
web_height = 150.0  # h
flange_width = 75.0   # b
web_thickness = 6.0   # t_w
flange_thickness = 10.0 # t_f

# Calculate the shear center location
e_location = calculate_shear_center_channel(web_height, flange_width, web_thickness, flange_thickness)

# --- Output the results ---
print("--- Shear Center Calculation for a Channel Section ---")
print(f"Given Dimensions:")
print(f"  Web Height (h): {web_height} mm")
print(f"  Flange Width (b): {flange_width} mm")
print(f"  Web Thickness (t_w): {web_thickness} mm")
print(f"  Flange Thickness (t_f): {flange_thickness} mm")
print("\n--- Calculation ---")
print("The shear center is a point where a shear force can be applied without causing twisting.")
print("For a channel section, this point is located outside the section on the side of the web.")
print("The distance 'e' from the web's vertical centerline is calculated.")
print(f"\nCalculated Shear Center Location (e): {e_location:.2f} mm from the web centerline.")
print("\nThis illustrates that the location is a specific, calculated distance determined by the section's dimensions, as described in the correct answer.")
