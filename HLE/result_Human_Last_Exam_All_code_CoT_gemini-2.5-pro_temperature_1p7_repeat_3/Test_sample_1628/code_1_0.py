import math

def calculate_channel_shear_center(h, b1, b2):
    """
    Calculates the shear center eccentricity for a thin-walled asymmetric channel section.

    The shear center is a geometric property of a cross-section. For an asymmetric
    channel section, it is located outside the section. Applying a shear force
    through this point will cause bending without twisting. This function calculates
    the eccentricity 'e' of the shear center from the web's centerline.

    Args:
        h (float): The height of the web (distance between flange centerlines).
        b1 (float): The width of the top flange.
        b2 (float): The width of the bottom flange.

    Returns:
        float: The eccentricity 'e' of the shear center from the web centerline.
    """
    # Based on thin-walled beam theory, the formula for eccentricity 'e' is:
    # e = (b1^2 + b2^2) / ( (2/3)*h + 2*(b1 + b2) )
    
    numerator = b1**2 + b2**2
    denominator = (2.0/3.0) * h + 2 * (b1 + b2)
    
    if denominator == 0:
        return float('inf')
        
    e = numerator / denominator
    
    # --- Output the step-by-step calculation ---
    print(f"Calculating the shear center eccentricity 'e' for a channel with:")
    print(f"  Web height (h) = {h}")
    print(f"  Top flange width (b1) = {b1}")
    print(f"  Bottom flange width (b2) = {b2}\n")
    
    print("The formula for the shear center eccentricity 'e' from the web centerline is:")
    print("e = (b1^2 + b2^2) / ((2/3)*h + 2*(b1 + b2))\n")
    
    print("Substituting the given values into the equation:")
    # Using .2f for cleaner display in the formula printout
    print(f"e = ({b1}^2 + {b2}^2) / ((2/3)*{h} + 2*({b1} + {b2}))")
    print(f"e = ({b1**2} + {b2**2}) / ({ (2.0/3.0) * h :.2f} + {2 * (b1 + b2)})")
    print(f"e = {numerator} / {denominator:.2f}\n")
    
    print(f"Result:")
    print(f"The shear center is located at a distance e = {e:.4f} units from the web centerline.")
    
    return e

# --- Example Usage ---
if __name__ == "__main__":
    # Define the dimensions of the asymmetric channel section
    web_height = 150.0  # e.g., in mm
    top_flange_width = 100.0  # e.g., in mm
    bottom_flange_width = 50.0   # e.g., in mm
    
    calculate_channel_shear_center(web_height, top_flange_width, bottom_flange_width)
