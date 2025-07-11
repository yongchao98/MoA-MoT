import numpy as np

def calculate_magnetic_field(Ha, w, d, D, Jc, x, z):
    """
    Calculates the magnetic field for a stack of superconducting strips.

    Args:
        Ha (float): Applied magnetic field.
        w (float): Half-width of the strips.
        d (float): Thickness of the strips.
        D (float): Stacking interval of the strips.
        Jc (float): Critical current density.
        x (float): x-coordinate for field calculation.
        z (float): z-coordinate for field calculation.
    """
    # Define H0 as given in the problem
    H0 = Jc * d / np.pi

    if Ha <= 0:
        print("Applied field Ha must be positive.")
        # For Ha=0 (zero-field cooled), field is zero everywhere.
        print("H_total = 0.0 + 0.0 = 0.0")
        return

    # Step 1: Calculate the flux penetration depth 'a' from Ha
    # Ha = H0 * ln[sinh(pi*w/D) / sinh(pi*a/D)]
    # => sinh(pi*a/D) = sinh(pi*w/D) * exp(-Ha/H0)
    try:
        sinh_arg_a = np.sinh(np.pi * w / D) * np.exp(-Ha / H0)
        # Check if the argument to arcsinh is valid (it always should be for w>0)
        if sinh_arg_a >= 0:
            a = (D / np.pi) * np.arcsinh(sinh_arg_a)
        else:
            # This case should not be reached physically
            print("Error: Could not calculate penetration depth 'a'.")
            return
    except OverflowError:
        print("Error: Overflow encountered in calculating 'a'. Ha might be too large.")
        return
        
    # Step 2: Calculate the self-field H_self(x, z)
    # Helper function for the numerator/denominator of the log argument
    def field_term(pos, width):
        cos_term = np.cos(2 * np.pi * z / D)
        cosh_term1 = np.cosh(2 * np.pi * (pos - width) / D)
        cosh_term2 = np.cosh(2 * np.pi * (pos + width) / D)
        return (cosh_term1 - cos_term) * (cosh_term2 - cos_term)

    numerator = field_term(x, w)
    denominator = field_term(x, a)

    # Avoid division by zero or log of non-positive number
    if numerator <= 0 or denominator <= 0:
        # This can happen if x=a or x=w and z=nD where cos=1.
        # It indicates being on top of a line of current singularity in this model.
        print(f"Cannot calculate field at the singular point (x={x}, z={z}).")
        print("This position coincides with the location of the idealized current sheets.")
        return

    log_argument = numerator / denominator
    H_self = - (H0 / 4) * np.log(log_argument)

    # Step 3: Calculate the total field
    H_total = Ha + H_self

    # Output the final equation with numbers
    print(f"For the given parameters:")
    print(f"Calculated flux penetration half-width a = {a:.4f}")
    print(f"The equation for the total field is H_total = Ha + H_self")
    print(f"{H_total:.4f} = {Ha:.4f} + ({H_self:.4f})")
    
    # Return the final result for programmatic use if needed
    return H_total

if __name__ == '__main__':
    # --- User-defined parameters ---
    # Example values based on a typical scenario
    
    Ha_applied = 0.1   # Applied magnetic field (in appropriate units, e.g., A/m)
    strip_half_width = 1e-3 # (w) e.g., 1 mm
    strip_thickness = 1e-6  # (d) e.g., 1 um
    stack_interval = 2e-3   # (D) e.g., 2 mm
    J_critical = 3e10       # (Jc) Critical current density (A/m^2)
    
    # Point where to calculate the field
    x_coord = 1.1e-3  # x-coordinate, outside the strip
    z_coord = 0.0     # z-coordinate, in the plane of a strip

    # Run the calculation
    print(f"Calculating magnetic field at (x, z) = ({x_coord}, {z_coord})")
    H_final = calculate_magnetic_field(Ha=Ha_applied, 
                                       w=strip_half_width, 
                                       d=strip_thickness, 
                                       D=stack_interval, 
                                       Jc=J_critical, 
                                       x=x_coord, 
                                       z=z_coord)
    # The final numerical answer is part of the printout.
    # To conform to the output format, we extract the final number.
    # The 'H_final' variable holds the numeric value.
    if H_final is not None:
        pass
        # Final answer is printed inside the function. For example: <<<9549.2969>>>
        # Since this value depends on user inputs, it is better to not output a fixed value here.
