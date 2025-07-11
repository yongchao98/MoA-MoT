import sympy
from sympy import init_printing

def solve_magnetostatics_problem():
    """
    This function presents the solution to the magnetostatics problem
    by symbolically representing and printing the final derived magnetic field H
    in both specified regions.
    """
    # Use pretty printing for the output
    init_printing(use_unicode=True)

    # Define symbolic variables
    M0, Rp, R, r, theta = sympy.symbols('M_0 R_p R r theta')
    
    # Define coordinate system basis vectors for clarity in printing
    ir = sympy.Symbol('î_r')
    itheta = sympy.Symbol('î_θ')

    print("This script displays the derived magnetic field H in the two regions.")
    print("--------------------------------------------------------------------\n")

    # --- Region 1: Inside the shield (0 < r < Rp) ---
    print("In the region 0 < r < R_p:")
    
    # The derived coefficient for the H-field inside the sphere.
    # The field is uniform in this region.
    coeff_H1 = M0 * (2*Rp**3 + R**3) / (3*R**3)
    
    # The vector part of the H-field in this region.
    # The tuple represents (r-component, theta-component)
    vector_part_1_str = "(-cos(θ) î_r + sin(θ) î_θ)"

    print("H(r, θ) = C * " + vector_part_1_str)
    print("where the coefficient C is:")
    sympy.pprint(coeff_H1, use_unicode=True)
    print("\nPutting it all together for Region 1:")
    print("H = (", coeff_H1, ") * ", vector_part_1_str)
    print("-" * 50)

    # --- Region 2: Between shield and conductor (Rp < r < R) ---
    print("\nIn the region R_p < r < R:")
    
    # Radial component of H-field (Hr)
    H2_r_coeff_part = (sympy.sympify(2)/3) * M0
    H2_r_radial_part = (Rp/R)**3 - (Rp/r)**3
    H2_r = -H2_r_coeff_part * H2_r_radial_part * sympy.cos(theta)
    
    # Theta component of H-field (Hθ)
    H2_theta_coeff_part = M0 / 3
    H2_theta_radial_part = 2*(Rp/R)**3 + (Rp/r)**3
    H2_theta = H2_theta_coeff_part * H2_theta_radial_part * sympy.sin(theta)
    
    print("The radial component Hr is:")
    print("Hr = (", -H2_r_coeff_part * sympy.cos(theta), ") * [", H2_r_radial_part, "]")
    print("\nThe tangential component Hθ is:")
    print("Hθ = (", H2_theta_coeff_part * sympy.sin(theta), ") * [", H2_theta_radial_part, "]")
    print("-" * 50)
    
    print("\nThese derived expressions match the equations given in answer choice B.")

if __name__ == '__main__':
    solve_magnetostatics_problem()
