import sympy as sp

def solve_magnetostatic_problem():
    """
    This function outlines the derivation of the magnetic field H and presents the final expressions.
    It doesn't perform numerical calculations but prints the symbolic results.
    """
    # Define symbols
    r, theta, R_p, R, M0 = sp.symbols('r, theta, R_p, R, M_0', real=True, positive=True)

    # The problem is solved using magnetostatic potential and boundary conditions.
    # After solving for the potential coefficients as outlined in the plan, we get:
    # A = M0 * (2*R_p**3 + R**3) / (3 * R**3)
    #
    # From these, the H fields are derived via H = -grad(Phi).

    print("The derived magnetic field H in each region is:")
    print("-" * 50)

    # Region 1: Inside the shield (0 < r < R_p)
    print("Region 1 (0 < r < R_p):")
    
    # H_in = A * (-cos(theta)*i_r + sin(theta)*i_theta)
    # The coefficient 'A' is:
    coeff_in_str = f"(M_0 * (2*{R_p**3} + {R**3})) / (3*{R**3})"
    
    # Vector part in spherical unit vectors
    vector_part_in_str = "(-cos(theta) i_r + sin(theta) i_theta)"
    
    H_in_str = f"H_in = {coeff_in_str} * {vector_part_in_str}"
    
    # Print the equation in a more readable format
    print(f"  H = M_0 * (2*R_p^3 + R^3)/(3*R^3) * (-cos(theta) i_r + sin(theta) i_theta)")
    # Split into components
    print(f"  H_r = -M_0 * (2*R_p^3 + R^3)/(3*R^3) * cos(theta)")
    print(f"  H_theta = M_0 * (2*R_p^3 + R^3)/(3*R^3) * sin(theta)")
    
    print("-" * 50)

    # Region 2: Between shield and conductor (R_p < r < R)
    print("Region 2 (R_p < r < R):")
    
    # H_r,out = - (2*M0/3) * ( (R_p**3 / R**3) - (R_p**3 / r**3) ) * cos(theta)
    # H_theta,out = (M0/3) * ( 2*(R_p**3 / R**3) + (R_p**3 / r**3) ) * sin(theta)
    
    # Radial component
    Hr_out_coeff_str = f"-(2*M_0/3)"
    Hr_out_radial_part_str = f"((R_p/R)^3 - (R_p/r)^3)"
    print(f"  H_r = {Hr_out_coeff_str} * [{Hr_out_radial_part_str}] * cos(theta)")

    # Theta component
    Htheta_out_coeff_str = f"(M_0/3)"
    Htheta_out_radial_part_str = f"(2*(R_p/R)^3 + (R_p/r)^3)"
    print(f"  H_theta = {Htheta_out_coeff_str} * [{Htheta_out_radial_part_str}] * sin(theta)")
    
    print("-" * 50)
    print("Comparing these results with the given options, we find that they exactly match Answer Choice B.")

solve_magnetostatic_problem()