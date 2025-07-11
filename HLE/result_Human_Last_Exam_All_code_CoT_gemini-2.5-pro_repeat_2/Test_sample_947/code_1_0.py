import math

def generate_field_expression():
    """
    This function derives and prints the expression for the magnetic field
    outside an infinite stack of superconducting strips.
    """

    # --- Symbolic representation of the equation's components ---

    # Main variables for the final equation
    # H_z(x,z) is the z-component of the magnetic field at position (x,z)
    field = "H_z(x,z)"
    # H_a is the externally applied magnetic field
    applied_field = "H_a"
    # H_0 is the characteristic field scale given by H0 = Jc * d / pi
    char_field = "H_0"

    # Physical parameters of the system
    # 2w is the width of the superconducting strips
    strip_half_width = "w"
    # D is the period of the stack (distance between strips)
    stack_period = "D"
    # a is the position of the flux front, determined by H_a
    flux_front = "a"
    
    # Position coordinates
    x_coord = "|x|"
    z_coord = "z"

    # Summation index
    index = "m"

    # --- Constructing the final expression as a string ---
    
    # This is the general term inside the infinite sum
    # Note: "2*pi*m" represents the m-th wavevector component k_m = 2*pi*m/D
    # exp(-2*pi*m*|x|/D) describes the exponential decay of the field away from the stack
    # [cosh(...) - cosh(...)] term relates to the geometry of the current distribution
    # cos(2*pi*m*z/D) describes the periodic variation of the field between the strips
    
    general_term = (
        f"(1/{index}) * exp(-2*pi*{index}*{x_coord}/{stack_period}) * "
        f"[cosh(2*pi*{index}*{strip_half_width}/{stack_period}) - cosh(2*pi*{index}*{flux_front}/{stack_period})] * "
        f"cos(2*pi*{index}*{z_coord}/{stack_period})"
    )

    # The full expression is the applied field minus the sum of contributions from all modes
    final_equation = (
        f"{field} = {applied_field} - {char_field} * Summation_{index}=1_to_inf [ {general_term} ]"
    )

    # --- Print the results ---
    
    print("Based on the Fourier analysis of the periodic current distribution, the expression for the magnetic field H_z is:")
    print("-" * 80)
    print(final_equation)
    print("-" * 80)
    print("\nWhere the terms are defined as:")
    print(f" H_z(x,z): The magnetic field's z-component at coordinates (x, z).")
    print(f" H_a:       The uniform applied magnetic field in the z-direction.")
    print(f" H_0:       The characteristic field, defined as H_0 = Jc*d/pi.")
    print(f" Jc:        The critical current density of the superconductor.")
    print(f" d:         The thickness of each superconducting strip.")
    print(f" pi:        The mathematical constant pi ({math.pi:.4f}...).")
    print(f" w:         Half the width of each strip (full width is 2w).")
    print(f" D:         The stacking interval between the strips.")
    print(f" a:         The parameter defining the flux-penetrated region (a < |x| < w). 'a' depends on H_a.")
    print(f" |x|, z:    The spatial coordinates. The expression is valid for |x| > w.")
    print(f" m:         The summation index for the Fourier modes (positive integers).")


# Execute the function to print the solution
generate_field_expression()

# The final expression in a more readable format is:
# H_z(x,z) = H_a - H_0 * Σ_{m=1 to ∞} (1/m) * exp(-2πm|x|/D) * [cosh(2πmw/D) - cosh(2πma/D)] * cos(2πmz/D)
final_answer = "H_z(x,z) = H_a - H_0 * Summation_{m=1_to_inf} [ (1/m) * exp(-2*pi*m*|x|/D) * [cosh(2*pi*m*w/D) - cosh(2*pi*m*a/D)] * cos(2*pi*m*z/D) ]"
print(f"\n<<<H_z(x,z) = H_a - H_0 * Summation_m=1_to_inf [ (1/m) * exp(-2*pi*m*|x|/D) * [cosh(2*pi*m*w/D) - cosh(2*pi*m*a/D)] * cos(2*pi*m*z/D) ]>>>")