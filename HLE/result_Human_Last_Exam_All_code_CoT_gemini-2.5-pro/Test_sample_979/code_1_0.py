def solve_magnetic_field():
    """
    This function prints the derived expressions for the magnetic field H
    inside and outside the spherical shell, based on the provided physics problem.
    """
    
    # Define symbolic constants from the problem description
    # K0: current amplitude
    # R: radius of the sphere
    # mu: magnetic permeability inside the sphere
    # mu0: magnetic permeability of free space (outside)
    # r: radial distance from the origin
    # theta: polar angle
    # H_in_str: string representing H-field inside the sphere
    # H_out_str: string representing H-field outside the sphere
    
    # Based on solving Laplace's equation with the given boundary conditions,
    # the magnetic field H is determined.
    
    # Inside the sphere (0 < r < R)
    # The field is uniform and directed along the z-axis.
    # The expression is derived as (2 * mu0 * K0) / (2 * mu0 + mu) * z_hat
    # which can be rewritten to match the format of the answer choices.
    H_in_str = "(2 * mu0 / mu) * K0 / (1 + (2 * mu0 / mu)) * z_hat"

    # Outside the sphere (r > R)
    # The field has the form of a magnetic dipole.
    # The expression is derived as (mu * R^3 * K0) / (2 * mu0 + mu) * (1/r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)
    # which can also be rewritten to match the format of the answer choices.
    H_out_str = "K0 / (1 + (2 * mu0 / mu)) * (R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)"

    print("The solution for the magnetic field H(r, theta) is:")
    print("-" * 50)
    
    print("For the region inside the sphere (0 < r < R):")
    print(f"H_in = {H_in_str}\n")
    
    print("For the region outside the sphere (r > R):")
    print(f"H_out = {H_out_str}\n")
    
    print("These results correspond to Answer Choice E.")
    
    # The problem requires printing the final equation with each number.
    # Since the result is symbolic, we print the components of the formulas.
    print("\n--- Final Equation from Choice E ---")
    print("H(r, theta) = ")
    print("{")
    # Inside field part
    print("  ( (2 * mu0) / mu ) * ( K0 / (1 + (2 * mu0 / mu)) ) * z_hat,  for 0 < r < R")
    # Outside field part
    print("  ( K0 / (1 + (2 * mu0 / mu)) ) * (R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat), for R < r < infinity")
    print("}")


solve_magnetic_field()