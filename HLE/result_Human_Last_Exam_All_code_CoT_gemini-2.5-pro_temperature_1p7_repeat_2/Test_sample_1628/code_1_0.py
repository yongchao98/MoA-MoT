import math

def calculate_channel_shear_center(depth, flange_width, web_thickness, flange_thickness):
    """
    Calculates the shear center location for a standard (singly symmetric) channel section.

    The shear center is the point through which an applied shear force
    produces no twisting moment. For a channel, it is located outside the
    section, on the opposite side of the centroid from the web.

    Args:
        depth (float): The total depth of the channel (d).
        flange_width (float): The width of the flanges (b_f).
        web_thickness (float): The thickness of the web (t_w).
        flange_thickness (float): The thickness of the flanges (t_f).
        
    Note: This calculation uses a common simplified formula for a thick-walled section and assumes sharp corners. 
    Results may differ slightly from tabulated values which account for corner radii.
    """
    print("--- Input Properties ---")
    print(f"Overall Depth (d): {depth}")
    print(f"Flange Width (bf): {flange_width}")
    print(f"Web Thickness (tw): {web_thickness}")
    print(f"Flange Thickness (tf): {flange_thickness}\n")

    # Step 1: Calculate the moment of inertia (I_x) about the horizontal centroidal axis
    # The centroidal axis is the axis of symmetry.
    
    # Clear height of the web
    h_clear = depth - 2 * flange_thickness
    # Distance between centroids of the two flanges
    h_centerline = depth - flange_thickness

    # Moment of inertia of the web about its own centroid + parallel axis theorem
    I_web = (web_thickness * h_clear**3) / 12
    
    # Moment of inertia of the two flanges about the section's centroidal axis
    # I_flange_own_axis + Area * distance^2
    I_flanges = 2 * ( (flange_width * flange_thickness**3 / 12) + \
                      (flange_width * flange_thickness) * (h_centerline / 2)**2 )
    
    I_x = I_web + I_flanges
    
    print("--- Intermediate Calculations ---")
    print(f"Distance between flange centerlines (h_centerline): {h_centerline:.4f}")
    print(f"Moment of Inertia about x-axis (Ix): {I_x:.4f}\n")
    
    # Step 2: Calculate the shear center location 'e' from the web centerline
    # The formula is derived from balancing the external shear force moment V*e with
    # the internal torque created by the shear flow in the flanges.
    # Formula: e = (tf * bf^2 * h_centerline^2) / (4 * Ix)
    
    e = (flange_thickness * flange_width**2 * h_centerline**2) / (4 * I_x)

    print("--- Final Calculation ---")
    print("Shear center distance 'e' from web centerline is calculated as:")
    print("e = (tf * bf^2 * h_centerline^2) / (4 * Ix)\n")
    
    # The problem asks to output each number in the final equation
    print("Substituting the values:")
    # Numerator calculation
    num = flange_thickness * flange_width**2 * h_centerline**2
    # Denominator calculation
    den = 4 * I_x
    print(f"e = ({flange_thickness:.4f} * {flange_width:.4f}^2 * {h_centerline:.4f}^2) / (4 * {I_x:.4f})")
    print(f"e = {num:.4f} / {den:.4f}\n")
    
    print("--- Result ---")
    print(f"The shear center is located at a distance 'e' = {e:.4f} from the web centerline.")

    # Often, the location is given from the back of the web (e0)
    e0 = e + web_thickness / 2
    print(f"The shear center is located at a distance 'e0' = {e0:.4f} from the back of the web.")

# Example using a standard AISC C10x15.3 channel section
# All units are in inches.
d_in = 10.0      # depth
bf_in = 2.60     # flange width
tw_in = 0.24     # web thickness
tf_in = 0.436    # flange thickness

calculate_channel_shear_center(d_in, bf_in, tw_in, tf_in)